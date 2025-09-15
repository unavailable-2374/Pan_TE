# RECON Modular Architecture

This directory contains the modular components of the advanced RECON pipeline, designed to improve maintainability, testability, and extensibility.

## Architecture Overview

The original monolithic `run_RECON_advanced` script (5700+ lines) has been refactored into focused, single-responsibility modules:

```
lib/RECON/
├── Config.pm         # Configuration management and argument parsing
├── Logger.pm         # Consistent logging with timing and resource tracking
├── Utils.pm          # Common utility functions (sampling, masking, BLAST)
├── Core.pm           # Core RECON algorithm implementation
├── MaskedTrack.pm    # High-quality consensus from masked genomes
├── SamplingTrack.pm  # Progressive sampling with adaptive stopping
└── README.md         # This documentation
```

## Module Responsibilities

### RECON::Config
- Command-line argument parsing and validation
- Input file discovery (genome, BED files)
- Resource allocation and memory management
- Processing mode determination (single vs dual track)
- Pipeline configuration based on genome size

### RECON::Logger
- Centralized logging with consistent formatting
- Pipeline timing and progress tracking
- Resource usage monitoring
- Step-based execution logging

### RECON::Utils
- Genome processing utilities (sampling, masking, indexing)
- BLAST database creation and management
- File format conversions and validations
- Size calculations and genome classification
- Tool availability checking

### RECON::Core
- Core RECON algorithm implementation
- Element definition, refinement, and family clustering
- K-parameter optimization
- RMBlastN self-alignment orchestration
- MSP (Maximal Scoring Pairs) collection

### RECON::MaskedTrack
- High-quality consensus building from masked genomes
- Unmasked region extraction for large genomes
- Full genome processing for small genomes
- Consensus sequence validation and quality control

### RECON::SamplingTrack
- Progressive multi-round sampling strategy
- Adaptive stopping criteria based on discovery metrics
- Cumulative masking with previous consensus libraries
- Round-based directory organization
- Metrics analysis and convergence detection

## Benefits of Modularization

### 1. **Maintainability**
- Each module has a single, well-defined responsibility
- Code is easier to understand, debug, and modify
- Changes to one component don't affect others
- Clear interfaces between components

### 2. **Testability**
- Individual modules can be unit tested
- Mock dependencies for isolated testing
- Easier to identify and fix bugs
- Regression testing of specific components

### 3. **Extensibility**
- New sampling strategies can be added to SamplingTrack
- Different logging outputs can be implemented in Logger
- Alternative masking algorithms in Utils
- Plugin architecture for custom processing modes

### 4. **Code Reuse**
- Utilities can be used by other Pan_TE components
- Core RECON functions available for standalone use
- Logger can be shared across the entire pipeline
- Configuration management reusable for other tools

## Usage Examples

### Basic Usage
```bash
# Same interface as before
bin/run_RECON_advanced 80 genome.fa 2500000000
```

### Programmatic Usage
```perl
use RECON::SamplingTrack;
use RECON::Logger;

# Run just the sampling track
run_sampling_track_independent($genome, \@bed_files, $consensi, 
                              $output_dir, $threads, $cpu_threads, $size);

# Use individual utilities
log_message("INFO", "Custom processing step", "details=example");
```

### Custom Configurations
```perl
use RECON::Config;

my $config = parse_arguments();
$config->{sample_sizes} = [20, 60, 180]; # Custom sampling sizes
```

## Migration Guide

### For Users
No changes required - the command-line interface remains identical:
```bash
# Old and new versions use the same syntax
bin/run_RECON_advanced [options] threads genome_file genome_size
```

### For Developers
The original monolithic version is preserved as:
```bash
bin/run_RECON_advanced.monolithic.backup
```

Key changes for developers:
1. Import required modules instead of defining functions
2. Use centralized logging through RECON::Logger
3. Access utilities through RECON::Utils
4. Configuration through RECON::Config

## Performance Considerations

The modular architecture maintains the same performance characteristics:
- **No overhead**: Perl's `use` statements are compile-time
- **Same algorithms**: Core logic unchanged, just better organized
- **Memory efficiency**: Modules loaded on-demand
- **Parallel processing**: Thread management preserved in each track

## Future Enhancements

The modular structure enables several planned improvements:

1. **Enhanced Sampling Strategies**
   - Coverage-based sampling
   - Phylogenetic-aware sampling
   - Graph genome sampling

2. **Advanced Stopping Criteria**
   - Machine learning-based convergence detection
   - Resource-aware adaptive stopping
   - Quality-based early termination

3. **Alternative Consensus Building**
   - Multiple sequence alignment improvements
   - Consensus validation pipelines
   - Quality scoring systems

4. **Integration Capabilities**
   - REST API for remote processing
   - Workflow management integration
   - Cloud computing adaptations

## Development Guidelines

When contributing to the modular RECON pipeline:

1. **Follow single responsibility principle**
2. **Use consistent error handling patterns**
3. **Add comprehensive logging for debugging**
4. **Write unit tests for new functions**
5. **Update documentation for interface changes**
6. **Maintain backward compatibility where possible**

## Dependencies

Each module specifies its own dependencies:
- **Config**: Getopt::Long, File::Spec
- **Logger**: Time::HiRes  
- **Utils**: File::Path, Cwd
- **Core**: (depends on Utils, Logger)
- **Tracks**: (depend on Utils, Core, Logger)

## Support

For questions about the modular architecture:
1. Check this README and module POD documentation
2. Review the examples in each module
3. Compare with the backup monolithic version
4. Submit issues to the Pan_TE repository

---

This modular architecture represents a significant improvement in code organization while maintaining full compatibility with existing workflows and performance expectations.