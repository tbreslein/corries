use super::DataValue;

pub struct Output<'a> {
    /// Index offset for printing spatially resolved data, useful for skipping ghost cells
    mesh_offset: usize,

    /// Whether this is this `Output`'s first output in this run
    first_output: bool,

    /// Whether this `Output` should write metadata into its `Stream`
    should_write_metadata: bool,

    /// Precision of floating point numbers
    precision: usize,

    /// Matrix of raw, unformatted data, still retaining their original data types
    data_matrix: Vec<&'a DataValue>,

    /// Intermediate string representation of the elements of the `data_matrix`
    string_matrix: Vec<Vec<String>>,

    /// Vector of strings ready to be written to a stream
    stream_strings: Vec<String>,

    /// Writes arbitrary data from `corries` objects to `string_matrix`
    string_conversion_mode: ToStringConversion,

    /// Formats the `string_matrix` into `stream_strings` which can be written to a `Stream`
    formatter: Formatter,

    /// Handles writing into a stream
    stream: Stream,
}

enum ToStringConversion {
    Scalar,
    Vector,
}

enum Formatter {
    CSV,
    TSV,
}

enum Stream {
    Stdout,
    File,
}
