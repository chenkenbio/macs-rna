//! Per-subcommand argument structs and entry points.
//!
//! Each module mirrors one MACS3 subcommand's argparse surface (from
//! `references/.../bin/macs3` + `Commands/<name>_cmd.py` +
//! `Utilities/OptValidator.py`). Field names use the MACS `dest` names so the
//! mapping is one-to-one. Bodies are `run(&Args) -> Result<()>` and currently
//! `todo!()`; the Phase 1 agents fill them.
//!
//! Ownership (frozen `mod.rs` after Phase 0): A owns `pileup`; B owns `bdgopt`,
//! `bdgpeakcall`, `bdgbroadcall`; C owns `bdgcmp`, `bdgdiff`; E owns `predictd`;
//! `callpeak` is Phase 2.

pub mod bdgbroadcall;
pub mod bdgcmp;
pub mod bdgdiff;
pub mod bdgopt;
pub mod bdgpeakcall;
pub mod callpeak;
pub mod pileup;
pub mod predictd;
