package RECON::Logger;

use strict;
use warnings;
use Exporter 'import';
use Time::HiRes qw(time);

our @EXPORT = qw(log_message log_step_start log_step_end log_progress log_resource_usage init_logger);

our $PIPELINE_START_TIME = time();
our $CURRENT_STEP_START = 0;
our $LOG_FILE_HANDLE = undef;
our $LOG_FILE_PATH = "recon_advanced.log";
our $LOG_FILE_INITIALIZED = 0;

sub log_message {
    my ($level, $message, $details) = @_;
    
    my $elapsed = time() - $PIPELINE_START_TIME;
    my $step_elapsed = $CURRENT_STEP_START > 0 ? time() - $CURRENT_STEP_START : 0;
    my $timestamp = sprintf("[%s] [%s] [+%.1fs] [step:%.1fs] %s", 
                           scalar(localtime), $level, $elapsed, $step_elapsed, $message);
    
    # Initialize log file handle if not already done
    if (!defined $LOG_FILE_HANDLE) {
        # Use overwrite mode on first initialization, append mode thereafter
        my $mode = $LOG_FILE_INITIALIZED ? '>>' : '>';
        open($LOG_FILE_HANDLE, $mode, $LOG_FILE_PATH) or die "Cannot open log file $LOG_FILE_PATH: $!\n";
        $LOG_FILE_HANDLE->autoflush(1);  # Ensure immediate flushing
        $LOG_FILE_INITIALIZED = 1;
    }
    
    # Output to both console and log file
    print "$timestamp\n";
    print $LOG_FILE_HANDLE "$timestamp\n";
    
    if ($details) {
        print "    Details: $details\n";
        print $LOG_FILE_HANDLE "    Details: $details\n";
    }
}

sub log_step_start {
    my ($step_name, $description) = @_;
    $description ||= "";
    
    $CURRENT_STEP_START = time();
    log_message("INFO", "=== STARTING: $step_name ===", $description);
}

sub log_step_end {
    my ($step_name, $success, $description) = @_;
    $description ||= "";
    
    my $step_duration = $CURRENT_STEP_START > 0 ? time() - $CURRENT_STEP_START : 0;
    my $status = $success ? "COMPLETED" : "FAILED";
    
    log_message("INFO", "=== $status: $step_name (${step_duration}s) ===", $description);
    $CURRENT_STEP_START = 0;
}

sub log_progress {
    my ($current, $total, $item_name) = @_;
    $item_name ||= "items";
    
    my $percent = $total > 0 ? int(($current / $total) * 100) : 0;
    log_message("INFO", "Progress: $current/$total $item_name", "${percent}% complete");
}

sub log_resource_usage {
    my ($checkpoint) = @_;
    $checkpoint ||= "checkpoint";
    
    # Get memory usage
    my $mem_info = `cat /proc/self/status | grep -E '^(VmPeak|VmSize|VmRSS)'`;
    chomp $mem_info;
    $mem_info =~ s/\n/, /g;
    
    # Get CPU load
    my $load_avg = `cat /proc/loadavg`;
    chomp $load_avg;
    my ($load1) = split /\s+/, $load_avg;
    
    log_message("DEBUG", "Resource usage at $checkpoint", "memory=($mem_info), load=$load1");
}

sub init_logger {
    # Reset logger state for new pipeline run
    $PIPELINE_START_TIME = time();
    $CURRENT_STEP_START = 0;
    $LOG_FILE_INITIALIZED = 0;
    
    # Close existing log file handle if open
    if (defined $LOG_FILE_HANDLE) {
        close($LOG_FILE_HANDLE);
        $LOG_FILE_HANDLE = undef;
    }
    
    # Create fresh log file (overwrite any existing)
    open(my $temp_fh, '>', $LOG_FILE_PATH) or die "Cannot create log file $LOG_FILE_PATH: $!\n";
    close($temp_fh);
}

1;

__END__

=head1 NAME

RECON::Logger - Logging utilities for RECON pipeline

=head1 SYNOPSIS

    use RECON::Logger;
    
    log_message("INFO", "Pipeline started", "threads=80");
    log_step_start("TRACK 1", "Processing masked genome");
    log_step_end("TRACK 1", 1, "Completed successfully");

=head1 DESCRIPTION

Provides consistent logging functionality for the RECON pipeline with timing,
progress tracking, and resource usage monitoring.

=cut