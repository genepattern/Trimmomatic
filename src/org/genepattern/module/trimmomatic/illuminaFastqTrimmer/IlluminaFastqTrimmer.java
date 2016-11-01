package org.genepattern.module.trimmomatic.illuminaFastqTrimmer;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import joptsimple.OptionParser;
import joptsimple.OptionSet;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.BooleanUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.usadellab.trimmomatic.Trimmomatic;

public class IlluminaFastqTrimmer {
    private static final String INPUT_FILE_1 = "input.file.1";
    private static final String INPUT_FILE_2 = "input.file.2";
    private static final String OUTPUT_FILENAME_BASE = "output.filename.base";
    private static final String ADAPTER_CLIP_SEQUENCE_FILE = "adapter.clip.sequence.file";
    private static final String ADAPTER_CLIP_SEED_MISMATCHES = "adapter.clip.seed.mismatches";
    private static final String ADAPTER_CLIP_PALINDROME_CLIP_THRESHOLD = "adapter.clip.palindrome.clip.threshold";
    private static final String ADAPTER_CLIP_SIMPLE_CLIP_THRESHOLD = "adapter.clip.simple.clip.threshold";
    private static final String ADAPTER_CLIP_MIN_LENGTH = "adapter.clip.min.length";
    private static final String ADAPTER_CLIP_KEEP_BOTH_READS = "adapter.clip.keep.both.reads";
    private static final String TRIM_LEADING_QUALITY_THRESHOLD = "trim.leading.quality.threshold";
    private static final String TRIM_TRAILING_QUALITY_THRESHOLD = "trim.trailing.quality.threshold";
    private static final String MAX_INFO_TARGET_LENGTH = "max.info.target.length";
    private static final String MAX_INFO_STRICTNESS = "max.info.strictness";
    private static final String SLIDING_WINDOW_SIZE = "sliding.window.size";
    private static final String SLIDING_WINDOW_QUALITY_THRESHOLD = "sliding.window.quality.threshold";
    private static final String MIN_READ_LENGTH = "min.read.length";
    private static final String EXTRA_STEPS = "extra.steps";
    private static final String CREATE_TRIMLOG = "create.trimlog";
    private static final String PHRED_ENCODING = "phred.encoding";
    private static final String CONVERT_PHRED_SCORES = "convert.phred.scores";
    private static final String THREADS = "threads";

    private static final Pattern BY_WHITESPACE_PATT = Pattern.compile("\\s+");

    public static OptionSet parseOptions(String[] options) {
        OptionParser optionParser = new OptionParser() {
            {
                accepts(INPUT_FILE_1).withRequiredArg();
                accepts(INPUT_FILE_2).withRequiredArg();
                accepts(OUTPUT_FILENAME_BASE).withRequiredArg();
                accepts(ADAPTER_CLIP_SEQUENCE_FILE).withRequiredArg();
                accepts(ADAPTER_CLIP_SEED_MISMATCHES).withRequiredArg().ofType(Integer.class);
                accepts(ADAPTER_CLIP_PALINDROME_CLIP_THRESHOLD).withRequiredArg().ofType(Integer.class);
                accepts(ADAPTER_CLIP_SIMPLE_CLIP_THRESHOLD).withRequiredArg().ofType(Integer.class);
                accepts(ADAPTER_CLIP_MIN_LENGTH).withRequiredArg().ofType(Integer.class);
                accepts(ADAPTER_CLIP_KEEP_BOTH_READS).withRequiredArg().ofType(Boolean.class);
                accepts(TRIM_LEADING_QUALITY_THRESHOLD).withRequiredArg().ofType(Integer.class);
                accepts(TRIM_TRAILING_QUALITY_THRESHOLD).withRequiredArg().ofType(Integer.class);
                accepts(MAX_INFO_TARGET_LENGTH).withRequiredArg().ofType(Integer.class);
                accepts(MAX_INFO_STRICTNESS).withRequiredArg().ofType(Float.class);
                accepts(SLIDING_WINDOW_SIZE).withRequiredArg().ofType(Integer.class);
                accepts(SLIDING_WINDOW_QUALITY_THRESHOLD).withRequiredArg().ofType(Integer.class);
                accepts(MIN_READ_LENGTH).withRequiredArg().ofType(Integer.class);
                accepts(EXTRA_STEPS).withRequiredArg();
                accepts(CREATE_TRIMLOG).withRequiredArg().ofType(Boolean.class);
                accepts(PHRED_ENCODING).withRequiredArg();
                accepts(CONVERT_PHRED_SCORES).withRequiredArg();
                accepts(THREADS).withRequiredArg().ofType(Integer.class);
            }
        };
        return optionParser.parse(options);
    }

    public IlluminaFastqTrimmer() {
    }

    /**
     * Build a processing pipeline based on the given set of options.  This includes the
     * basic general options, the input and output files, and the processing operations.
     */
    public List<String> buildProcessingPipeline(OptionSet optionSet) throws TrimmingException {
        // This whole process is a little awkward.  The entry point into Trimmomatic is 
        // the main() method of a Trimmomatic instance.  That means that we need to build
        // a String array equivalent to the expected command line and pass that through.
        List<String> trimmomaticArgs = new ArrayList<String>();
        addBaseArgs(optionSet, trimmomaticArgs);

        // Get the length of the args in order to check whether the following add any steps
        int baseArgsLength = trimmomaticArgs.size();

        addConvertQualStep(optionSet, trimmomaticArgs);
        addIlluminaClipStep(optionSet, trimmomaticArgs);
        addIntegerOnlyStep(optionSet, TRIM_LEADING_QUALITY_THRESHOLD, "LEADING", trimmomaticArgs);
        addIntegerOnlyStep(optionSet, TRIM_TRAILING_QUALITY_THRESHOLD, "TRAILING", trimmomaticArgs);
        addMaxInfoStep(optionSet, trimmomaticArgs);
        addSlidingWindowStep(optionSet, trimmomaticArgs);
        addIntegerOnlyStep(optionSet, MIN_READ_LENGTH, "MINLEN", trimmomaticArgs);
        addExtraSteps(optionSet, trimmomaticArgs);

        // Now check to see whether any processing options were specified
        if (baseArgsLength == trimmomaticArgs.size()) {
            throw new TrimmingException("No processing steps were specified.  Please specify at least one operation.");
        }

        return trimmomaticArgs;
    }

    /**
     * Execute Trimmomatic with the given processing pipeline.  This handles error checking as well.
     */
    public void executeProcessingPipeline(List<String> trimmomaticArgs) throws TrimmingException {
        final String[] args = trimmomaticArgs.toArray(new String[trimmomaticArgs.size()]);
        printCommandLine(args);

        try {
            Trimmomatic.main(args);
        }
        catch (IOException ie) {
            throw new TrimmingException("A problem occurred during read trimming", ie);
        }
    }

    /**
     * Add the basic args and input/output files to the cmd-line based on the given OptionSet.  This
     * handles everything up to the actual processing steps.
     */
    private void addBaseArgs(OptionSet optionSet, List<String> trimmomaticArgs) throws TrimmingException {
        // Check the required parameters not governed by drop-downs 
        if (!optionSet.has(INPUT_FILE_1)) {
            throw new TrimmingException("Missing required parameter " + INPUT_FILE_1);
        }
        if (!optionSet.has(OUTPUT_FILENAME_BASE)) {
            throw new TrimmingException("Missing required parameter " + OUTPUT_FILENAME_BASE);
        }
        boolean isPairedEnded = optionSet.has(INPUT_FILE_2);
        trimmomaticArgs.add((isPairedEnded) ? "PE" : "SE");

        if (optionSet.has(THREADS)) {
            trimmomaticArgs.add("-threads");
            trimmomaticArgs.add(getNonnegativeIntegerValue(optionSet, THREADS));
        }

        String outputFilenameBase = optionSet.valueOf(OUTPUT_FILENAME_BASE).toString();
        if (optionSet.has(PHRED_ENCODING)) {
            String phredEncoding = optionSet.valueOf(PHRED_ENCODING).toString();
            if (StringUtils.equalsIgnoreCase(phredEncoding, "phred33")) {
                trimmomaticArgs.add("-phred33");
            }
            else if (StringUtils.equalsIgnoreCase(phredEncoding, "phred64")) {
                trimmomaticArgs.add("-phred64");
            }
        }

        if (getBoolean(optionSet, CREATE_TRIMLOG)) {
            trimmomaticArgs.add("-trimlog");
            trimmomaticArgs.add(outputFilenameBase + ".trimlog.txt");
        }

        String inputFile1 = getString(optionSet, INPUT_FILE_1);
        outputFilenameBase = cleanOutputFilenameBase(outputFilenameBase);
        String outputFilename = getOutputFilename(outputFilenameBase, inputFile1, isPairedEnded);

        if (isPairedEnded) {
            trimmomaticArgs.add("-baseout");
            trimmomaticArgs.add(outputFilename);

            trimmomaticArgs.add(inputFile1);
            trimmomaticArgs.add(getString(optionSet, INPUT_FILE_2));
        }
        else {
            trimmomaticArgs.add(inputFile1);
            trimmomaticArgs.add(outputFilename);
        }
    }

    // As the outputFilenameBase is <input.file.1_basename> by default, it can include various
    // undesirable filename-based artifacts due to the compression extension or the "_1" common
    // many of these fastqs.  Trim these off if possible.
    private String cleanOutputFilenameBase(String outputFilenameBase) {
        String cleaned = outputFilenameBase;
        cleaned = trimFilenameEndingIfNecessary(cleaned, ".fq");
        cleaned = trimFilenameEndingIfNecessary(cleaned, ".fastq");
        cleaned = trimFilenameEndingIfNecessary(cleaned, "_1");
        return cleaned;
    }

    private String trimFilenameEndingIfNecessary(String filename, String ending) {
        // Note that we check the length as we go to make sure that we aren't left with a blank name.
        int filenameLength = StringUtils.length(filename);
        int endingLength = StringUtils.length(ending);
        if (filenameLength > endingLength && StringUtils.endsWithIgnoreCase(filename, ending)) {
            return StringUtils.substring(filename, 0, filenameLength - endingLength);
        }
        return filename;
    }

    /**
     * Create an output filename suitable for use in either PE or SE mode.  
     * Trimmomatic will automatically compress the output based on the file extension.
     * Also, we want to match the .fq or .fastq usage in the original input file.
     * Note that for PE mode this is the string used by "-baseout" and not in and of
     * itself the complete filename.
     */
    private String getOutputFilename(String outputFilenameBase, String inputFile, boolean isPairedEnded) {
        String ext_compr = "";
        String ext = FilenameUtils.getExtension(inputFile);

        // Does this inputFile indicate compression?  If so, do likewise on the output
        if (StringUtils.endsWithIgnoreCase(inputFile, "gz") || StringUtils.endsWithIgnoreCase(inputFile, "bz2")) {
            // NOTE: Trimmomatic 0.32 has a major bug in compressing output with bz2; it throws an exception
            // during compression but swallows this and hangs indefinitely if run multithreaded (which we do).
            // If run single threaded it just fails.  Since this makes it unusable, we're always going to
            // provide compressed output as gz.
            ext_compr = ".gz";

            // Disabled for now because of the above issue
            // Save this as the "compression" extension.
            //ext_compr = "." + ext;

            // What's the extension when ext_compr is removed?
            String filenameWithNoExt = FilenameUtils.getBaseName(inputFile);
            ext = FilenameUtils.getExtension(filenameWithNoExt);
        }

        // If the input file has a recognized extension then match that for the output, otherwise use fq.
        if (!StringUtils.equalsIgnoreCase(ext, "fq") && !StringUtils.equalsIgnoreCase(ext, "fastq")) {
            ext = "fq";
        }

        if (isPairedEnded) {
            return outputFilenameBase + "." + ext + ext_compr;
        }

        // For single-ended input, make sure that the constructed filename is different from the inputFile
        return outputFilenameBase + "_trimmed." + ext + ext_compr;
    }

    /**
     * Print the equivalent command-line for debug/test purposes.
     */
    private void printCommandLine(String[] args) {
        try {
            File cmdlineFile = new File(".", "cmdline.log");
            FileWriter cmdlineWriter = new FileWriter(cmdlineFile);
            try {
                cmdlineWriter.write("Equivalent command line:");
                cmdlineWriter.write(IOUtils.LINE_SEPARATOR);
                cmdlineWriter.write("Trimmomatic");
                for (String arg : args) {
                    cmdlineWriter.write(" '");
                    cmdlineWriter.write(arg);
                    cmdlineWriter.write("'");
                }
                cmdlineWriter.write(IOUtils.LINE_SEPARATOR);
            }
            finally {
                if (cmdlineWriter != null) {
                    cmdlineWriter.close();
                }
            }
        }
        catch (IOException ie) {

        }
    }

    /**
     * Convenience method to get a simple String option value
     */
    private String getString(OptionSet optionSet, String optionName) {
        return optionSet.valueOf(optionName).toString();
    }

    /** 
     * Check that the given option value is a non-negative integer.  
     * Returns a String since we don't actually use the numeric value.
     *  
     **/
    private String getNonnegativeIntegerValue(OptionSet optionSet, String optionName) throws TrimmingException {
        if (!optionSet.has(optionName)) {
            return null;
        }

        Integer value = (Integer) optionSet.valueOf(optionName);
        if (value < 0) {
            throw new TrimmingException("Parameter '" + optionName + "' must be a non-negative integer.");
        }

        return value.toString();
    }

    /**
     * Build a step that requires only a single integer value.
     */
    private void addIntegerOnlyStep(OptionSet optionSet, String optionName, String stepCommand,
            List<String> trimmomaticArgs)
            throws TrimmingException {
        String value = getNonnegativeIntegerValue(optionSet, optionName);
        if (StringUtils.isBlank(value)) {
            return;
        }

        trimmomaticArgs.add(stepCommand + ":" + value);
    }

    /** 
     * Get the given option value as a boolean (as handled by jopt-simple).  
     **/
    private boolean getBoolean(OptionSet optionSet, String optionName) {
        return optionSet.has(optionName) && (Boolean) optionSet.valueOf(optionName);
    }

    /** 
     * Check that the given option value is a non-negative float.
     * Returns a String since we don't actually use the numeric value.
     *  
     **/
    private String getZeroToOneFloatValue(OptionSet optionSet, String optionName) throws TrimmingException {
        if (!optionSet.has(optionName)) {
            return null;
        }

        String token = optionSet.valueOf(optionName).toString();
        float value = NumberUtils.toFloat(token, -1);
        if (value < 0 || value > 1) {
            throw new TrimmingException("Parameter '" + optionName + "' must be a decimal value between 0 and 1.");
        }

        return token;
    }

    private void addConvertQualStep(OptionSet optionSet, List<String> trimmomaticArgs) {
        if (optionSet.has(CONVERT_PHRED_SCORES)) {
            String phredConversion = optionSet.valueOf(CONVERT_PHRED_SCORES).toString();
            if (StringUtils.equalsIgnoreCase(phredConversion, "to_phred33")) {
                trimmomaticArgs.add("TOPHRED33");
            }
            else if (StringUtils.equalsIgnoreCase(phredConversion, "to_phred64")) {
                trimmomaticArgs.add("TOPHRED64");
            }
        }
    }

    private void addIlluminaClipStep(OptionSet optionSet, List<String> trimmomaticArgs) throws TrimmingException {
        boolean hasSequenceFile = optionSet.has(ADAPTER_CLIP_SEQUENCE_FILE);
        boolean hasSeedMismatches = optionSet.has(ADAPTER_CLIP_SEED_MISMATCHES);
        boolean hasPalindromeClipThreshold = optionSet.has(ADAPTER_CLIP_PALINDROME_CLIP_THRESHOLD);
        boolean hasSimpleClipThreshold = optionSet.has(ADAPTER_CLIP_SIMPLE_CLIP_THRESHOLD);
        boolean hasMinLength = optionSet.has(ADAPTER_CLIP_MIN_LENGTH);
        boolean keepBothReads = getBoolean(optionSet, ADAPTER_CLIP_KEEP_BOTH_READS);
        boolean[] required = {
                hasSequenceFile, hasSeedMismatches, hasPalindromeClipThreshold, hasSimpleClipThreshold
        };

        // Check whether the user provided any/all of the parameters required for ILLUMINACLIP.  Also
        // check for the optional ADAPTER_CLIP_MIN_LENGTH and ADAPTER_CLIP_KEEP_BOTH_READS just to 
        // give a sensible error message if req. are missing.  For ADAPTER_CLIP_KEEP_BOTH_READS, we
        // bail if it is true (the default) but continue if it has been changed.
        if (!BooleanUtils.or(required) && !hasMinLength && keepBothReads) {
            // None of the options were provided.
            return;
        }
        if (!BooleanUtils.and(required)) {
            // Some options were provided, but some of the required are missing.
            throw new TrimmingException("For adapter clipping, you must provide adapter.clip.sequence.file, " +
                    "adapter.clip.seed.mismatches, adapter.clip.palindrome.clip.threshold, " +
                    "and adapter.clip.simple.clip.threshold.");
        }

        StringBuilder builder = new StringBuilder("ILLUMINACLIP");
        builder.append(":").append(optionSet.valueOf(ADAPTER_CLIP_SEQUENCE_FILE));
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, ADAPTER_CLIP_SEED_MISMATCHES));
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, ADAPTER_CLIP_PALINDROME_CLIP_THRESHOLD));
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, ADAPTER_CLIP_SIMPLE_CLIP_THRESHOLD));

        if (hasMinLength) {
            builder.append(":").append(getNonnegativeIntegerValue(optionSet, ADAPTER_CLIP_MIN_LENGTH));
            if (keepBothReads) {
                builder.append(":").append("true");
            }
        }
        else if (keepBothReads) {
            // ADAPTER_CLIP_MIN_LENGTH was not specified, but since we have to provide something we pass the default.
            // This is only necessary if we are also specifying keepBothReads
            builder.append(":").append("8");
            builder.append(":").append("true");
        }

        trimmomaticArgs.add(builder.toString());
    }

    private void addSlidingWindowStep(OptionSet optionSet, List<String> trimmomaticArgs) throws TrimmingException {
        boolean hasWindowSize = optionSet.has(SLIDING_WINDOW_SIZE);
        boolean hasQualityThreshold = optionSet.has(SLIDING_WINDOW_QUALITY_THRESHOLD);
        boolean[] required = { hasWindowSize, hasQualityThreshold };

        // Check whether the user provided any/all of the parameters required for SLIDINGWINDOW
        if (!BooleanUtils.or(required)) {
            // None of the required options were provided.
            return;
        }
        if (!BooleanUtils.and(required)) {
            // Only some but not all of the required options were provided.
            throw new TrimmingException(
                    "For the sliding window trim, you must provide sliding.window.quality.threshold and sliding.window.size");
        }

        StringBuilder builder = new StringBuilder("SLIDINGWINDOW");
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, SLIDING_WINDOW_SIZE));
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, SLIDING_WINDOW_QUALITY_THRESHOLD));
        trimmomaticArgs.add(builder.toString());
    }

    private void addMaxInfoStep(OptionSet optionSet, List<String> trimmomaticArgs) throws TrimmingException {
        boolean hasTargetLength = optionSet.has(MAX_INFO_TARGET_LENGTH);
        boolean hasStrictness = optionSet.has(MAX_INFO_STRICTNESS);
        boolean[] required = { hasTargetLength, hasStrictness };

        // Check whether the user provided any/all of the parameters required for MAXINFO
        if (!BooleanUtils.or(required)) {
            // None of the required options were provided.
            return;
        }
        if (!BooleanUtils.and(required)) {
            // Only some but not all of the required options were provided.
            throw new TrimmingException(
                    "For the max info trim, you must provide max.info.target.length and max.info.strictness");
        }

        StringBuilder builder = new StringBuilder("MAXINFO");
        builder.append(":").append(getNonnegativeIntegerValue(optionSet, MAX_INFO_TARGET_LENGTH));
        builder.append(":").append(getZeroToOneFloatValue(optionSet, MAX_INFO_STRICTNESS));
        trimmomaticArgs.add(builder.toString());
    }

    private void addExtraSteps(OptionSet optionSet, List<String> trimmomaticArgs) throws TrimmingException {
        if (!optionSet.has(EXTRA_STEPS)) {
            return;
        }
        String extraStepValue = getString(optionSet, EXTRA_STEPS);
        String[] extraSteps = BY_WHITESPACE_PATT.split(extraStepValue);
        if (extraSteps == null) {
            return;
        }

        // We add these to the args as they are and just let Trimmomatic deal with them.
        for (String extraStep : extraSteps) {
            if (StringUtils.isNotBlank(extraStep)) {
                trimmomaticArgs.add(extraStep);
            }
        }
    }

    // Set up a default UncaughtExceptionHandler to deal with any unchecked exceptions thrown by the Trimmomatic
    // worker threads.  Trimmomatic is not handling them in a way that allows GP to detect that something went wrong. 
    // These will be reported on System.err so that the run will be flagged as an error; see notes below for the
    // reason that we need to hang on to System.err.
    private static PrintStream systemErr = System.err;
    private static TrimmingUncaughtExceptionHandler uncaughtExceptionHandler = new TrimmingUncaughtExceptionHandler(
            systemErr);
    static {
        Thread.setDefaultUncaughtExceptionHandler(uncaughtExceptionHandler);
    }

    public static void main(String[] args) throws TrimmingException {
        OptionSet optionSet = parseOptions(args);
        IlluminaFastqTrimmer trimmer = new IlluminaFastqTrimmer();
        List<String> trimmomaticArgs = trimmer.buildProcessingPipeline(optionSet);
        
        // Trimmomatic prints messages from normal operations to stderr which of course causes GP
        // to treat every run as an error.  Here we redirect STDERR to STDOUT (internally) to 
        // work around the issue. Fortunately, Trimmomatic seems to throw Exceptions as 
        // well when it hits a problem; these are dealt with by the uncaughtExceptionHandler.
        // NOTE: Trimmomatic in certain cases calls System.exit(1) when it encounters an error.
        // This effectively stops any post-processing, so we can't count on anything happening after
        // this point. It would be better to spawn a new JVM process, wait for that to complete and
        // then do any post-processing.  Fortunately, at this point it seems that nothing of that
        // sort will be necessary.

        try {
            // Redirect STDERR to STDOUT for the Trimmomatic call.
            System.setErr(System.out);
            trimmer.executeProcessingPipeline(trimmomaticArgs);
        }
        finally {
            // Restore the normal System.err.  Note that this may fail in cases where Trimmomatic calls
            // System.exit(1) directly rather than using an Exception.
            System.setErr(systemErr);
        }
    }
}
