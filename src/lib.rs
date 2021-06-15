use std::fmt;
use std::fs::File;
use std::collections::HashMap;

extern crate bio;
use bio::io::fastq;
use bio::io::fastq::FastqRead;

extern crate chrono;
use chrono::{DateTime, FixedOffset};

extern crate permutation;


#[allow(non_snake_case)]
pub struct Dataframe
{
    seq_length   : Vec<usize>,
    mean_quality : Vec<f64>,
    kmer_start   : Vec<[u8;4]>,
    kmer_end     : Vec<[u8;4]>,
    ntc_A        : Vec<usize>,
    ntc_G        : Vec<usize>,
    ntc_T        : Vec<usize>,
    ntc_C        : Vec<usize>,
    ntc_U        : Vec<usize>,
    channel      : Vec<usize>,
    start_time   : Vec<DateTime<FixedOffset>>,
}

impl fmt::Display for Dataframe
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result
    {
        write!(f, "seq_length\tmean_quality\tkmers_start\tkmers_end\tnt_A\tnt_G\tnt_T\tnt_C\tnt_U\tchannels\tstart_times\n")?;
        for i in 0..self.seq_length.len()
        {
            write!(f, "{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                   self.seq_length[i],
                   self.mean_quality[i],
                   String::from_utf8(self.kmer_start[i].to_vec()).expect("kmer_start failed"),
                   String::from_utf8(self.kmer_end[i].to_vec()).expect("kmer_end failed"),
                   self.ntc_A[i],
                   self.ntc_G[i],
                   self.ntc_T[i],
                   self.ntc_C[i],
                   self.ntc_U[i],
                   self.channel[i],
                   self.start_time[i],
                   )?;
        }
        return Ok(());
    }
}

pub fn write_df(df : Dataframe)
{
    println!("{}", df);
}

use std::io;
#[allow(non_snake_case)]
pub fn fastq_df(path: String) -> io::Result<String>
{
    let mut seq_length   = Vec::new();
    let mut mean_quality = Vec::new();
    let mut kmer_start   = Vec::new();
    let mut kmer_end     = Vec::new();
    let mut ntc_A        = Vec::new();
    let mut ntc_G        = Vec::new();
    let mut ntc_T        = Vec::new();
    let mut ntc_C        = Vec::new();
    let mut ntc_U        = Vec::new();
    let mut channel      = Vec::new();
    let mut start_time   = Vec::new();

    let file = File::open(&path)?;
    let mut reader = fastq::Reader::new(file);
    let mut record = fastq::Record::new();

    loop
    {
        reader.read(&mut record).expect("read error");
        if record.is_empty()
        {
            break
        }

        seq_length.push(record.seq().len());

        let k = if record.seq().len() < 4 { record.seq().len() } else { 4 };

        let mut kmer : [u8;4] = [0,0,0,0];
        kmer[0..k].copy_from_slice(&record.seq()[0..k]);
        kmer_start.push(kmer);

        let mut kmer : [u8;4] = [0,0,0,0];
        kmer[0..k].copy_from_slice(&record.seq()[record.seq().len()-k..]);
        kmer_end.push(kmer);

        let sum = record.qual().iter().fold(0.0, |sum, x| sum + (*x - 33) as f64);
        mean_quality.push(sum / record.seq().len() as f64);

        let mut nts : [usize;256] = [0;256];
        for b in record.seq()
        {
            nts[*b as usize] += 1;
        }

        ntc_A.push(nts['A' as usize]);
        ntc_G.push(nts['G' as usize]);
        ntc_T.push(nts['T' as usize]);
        ntc_C.push(nts['C' as usize]);
        ntc_U.push(nts['U' as usize]);

        let info : HashMap<String, String> = record.desc().unwrap().split(' ')
            .map(|kv| kv.split('=').collect::<Vec<&str>>())
            .map(|vec| { (vec[0].to_string(), vec[1].to_string()) })
            .collect();

        channel.push(info["ch"].parse().unwrap());
        let st = DateTime::parse_from_rfc3339(&info["start_time"]);
        match st
        {
            Ok(time) => start_time.push(time),
            Err(e) => eprintln!("Error parsing time: {}", e)
        }
    }

    let order = permutation::sort(&start_time[..]);
    let df = Dataframe
    {
        seq_length :   order.apply_slice(&seq_length[..]),
        mean_quality:  order.apply_slice(&mean_quality[..]),
        kmer_start:    order.apply_slice(&kmer_start[..]),
        kmer_end:      order.apply_slice(&kmer_end[..]),
        ntc_A:         order.apply_slice(&ntc_A[..]),
        ntc_G:         order.apply_slice(&ntc_G[..]),
        ntc_T:         order.apply_slice(&ntc_T[..]),
        ntc_C:         order.apply_slice(&ntc_C[..]),
        ntc_U:         order.apply_slice(&ntc_U[..]),
        channel:       order.apply_slice(&channel[..]),
        start_time:    order.apply_slice(&start_time[..])
    };
    return Ok(df.to_string());
}


#[cfg(test)]
mod tests
{
    #[test]
    fn it_works()
    {
        assert_eq!(2 + 2, 4);
    }
}

