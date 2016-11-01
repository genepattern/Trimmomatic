package org.genepattern.module.trimmomatic.illuminaFastqTrimmer;

public class TrimmingException extends RuntimeException {

    public TrimmingException() {
        super();
    }

    public TrimmingException(String cause, Throwable message) {
        super(cause, message);
    }

    public TrimmingException(String cause) {
        super(cause);
    }

    public TrimmingException(Throwable cause) {
        super(cause);
    }
}
