use lrdf::fastq_df;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyfunction]
pub fn fastq_df_py(path: String) -> PyResult<String> 
{
    Ok(fastq_df(path).expect(""))
}

#[pymodule]
fn lrdf(_py: Python, m: &PyModule) -> PyResult<()>
{
    //m.add_wrapped(wrap_pyfunction!(fastq_df))?;
    m.add_function(wrap_pyfunction!(fastq_df_py, m)?)?;
    Ok(())
}

