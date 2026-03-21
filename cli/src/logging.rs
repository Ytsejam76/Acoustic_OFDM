// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;

use log::LevelFilter;
use simplelog::{
    ColorChoice, CombinedLogger, ConfigBuilder, SharedLogger, TermLogger, TerminalMode,
    WriteLogger,
};

pub fn init_logging(log_file: Option<&str>, verbose: bool) -> Result<(), Box<dyn Error>> {
    let level = if verbose {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    let cfg = ConfigBuilder::new().set_time_format_rfc3339().build();
    let mut loggers: Vec<Box<dyn SharedLogger>> = Vec::new();
    loggers.push(TermLogger::new(
        level,
        cfg.clone(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    ));
    if let Some(path) = log_file {
        let file = std::fs::File::create(path)?;
        loggers.push(WriteLogger::new(level, cfg, file));
    }
    let _ = CombinedLogger::init(loggers);
    Ok(())
}

#[macro_export]
macro_rules! info_line {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        println!("{}", msg);
        log::info!("{}", msg);
    }};
}

#[macro_export]
macro_rules! warn_line {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        eprintln!("{}", msg);
        log::warn!("{}", msg);
    }};
}

// vim: set ts=4 sw=4 et:
