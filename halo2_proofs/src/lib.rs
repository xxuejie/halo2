//! # halo2_proofs

#![no_std]
#![cfg_attr(docsrs, feature(doc_cfg))]
// The actual lints we want to disable.
#![allow(clippy::op_ref, clippy::many_single_char_names)]
#![deny(rustdoc::broken_intra_doc_links)]
#![deny(missing_debug_implementations)]
#![deny(missing_docs)]
#![deny(unsafe_code)]

extern crate alloc;

use hashbrown as collections;

pub mod arithmetic;
pub mod circuit;
pub use pasta_curves as pasta;
pub mod plonk;
pub mod poly;
pub mod transcript;

mod helpers;
mod io;
