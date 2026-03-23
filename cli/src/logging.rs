// Copyright (c) 2026 Elias S. G. Carotti

use std::error::Error;

use clap::ValueEnum;
use log::LevelFilter;
use simplelog::{
    ColorChoice, CombinedLogger, ConfigBuilder, SharedLogger, TermLogger, TerminalMode, WriteLogger,
};

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum LogLevelArg {
    Error,
    Warn,
    Info,
    Debug,
    Trace,
}

impl From<LogLevelArg> for LevelFilter {
    fn from(value: LogLevelArg) -> Self {
        match value {
            LogLevelArg::Error => LevelFilter::Error,
            LogLevelArg::Warn => LevelFilter::Warn,
            LogLevelArg::Info => LevelFilter::Info,
            LogLevelArg::Debug => LevelFilter::Debug,
            LogLevelArg::Trace => LevelFilter::Trace,
        }
    }
}

pub fn init_logging(log_file: Option<&str>, level: LevelFilter) -> Result<(), Box<dyn Error>> {
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
        log::info!("{msg}");
    }};
}

#[macro_export]
macro_rules! warn_line {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        log::warn!("{msg}");
    }};
}

#[macro_export]
macro_rules! debug_line {
    ($($arg:tt)*) => {{
        let msg = format!($($arg)*);
        log::debug!("{msg}");
    }};
}

// vim: set ts=4 sw=4 et:
