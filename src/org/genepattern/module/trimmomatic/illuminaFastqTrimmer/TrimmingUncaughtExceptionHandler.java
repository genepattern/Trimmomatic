package org.genepattern.module.trimmomatic.illuminaFastqTrimmer;

import java.io.PrintStream;
import java.lang.Thread.UncaughtExceptionHandler;

public class TrimmingUncaughtExceptionHandler implements UncaughtExceptionHandler {

    private PrintStream systemErr;

    public TrimmingUncaughtExceptionHandler(PrintStream systemErr) {
        this.systemErr = systemErr;
    }

    @Override
    public void uncaughtException(Thread thread, Throwable throwable) {
        // Print the stackTrace of any uncaught exception to the *real* STDERR.  Note that this was redirected 
        // to System.out in general - here we want to force it to go to the process's STDERR.
        throwable.printStackTrace(systemErr);
    }
}
