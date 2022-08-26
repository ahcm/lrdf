use structopt::StructOpt;
use lrdf::fastq_df;
use lrdf::write_df;

#[derive(Debug, StructOpt)]
#[structopt(name = "lrdf", about = "Generate a data frame from ONT long read Fastq files")]
struct Opt
{
    //#[structopt(parse(from_os_str))]
    input: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>>
{
    let opt = Opt::from_args();
    let df = fastq_df(opt.input)?;
    Ok(write_df(df))
}
