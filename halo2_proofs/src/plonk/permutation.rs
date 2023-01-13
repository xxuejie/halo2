use super::circuit::{Any, Column};
use crate::{arithmetic::CurveAffine, helpers::CurveRead, io};
use alloc::{vec, vec::Vec};

pub(crate) mod keygen;
pub(crate) mod verifier;

/// A permutation argument.
#[derive(Debug, Clone)]
pub(crate) struct Argument {
    /// A sequence of columns involved in the argument.
    columns: Vec<Column<Any>>,
}

impl Argument {
    pub(crate) fn new() -> Self {
        Argument { columns: vec![] }
    }

    /// Returns the minimum circuit degree required by the permutation argument.
    /// The argument may use larger degree gates depending on the actual
    /// circuit's degree and how many columns are involved in the permutation.
    pub(crate) fn required_degree(&self) -> usize {
        // degree 2:
        // l_0(X) * (1 - z(X)) = 0
        //
        // We will fit as many polynomials p_i(X) as possible
        // into the required degree of the circuit, so the
        // following will not affect the required degree of
        // this middleware.
        //
        // (1 - (l_last(X) + l_blind(X))) * (
        //   z(\omega X) \prod (p(X) + \beta s_i(X) + \gamma)
        // - z(X) \prod (p(X) + \delta^i \beta X + \gamma)
        // )
        //
        // On the first sets of columns, except the first
        // set, we will do
        //
        // l_0(X) * (z(X) - z'(\omega^(last) X)) = 0
        //
        // where z'(X) is the permutation for the previous set
        // of columns.
        //
        // On the final set of columns, we will do
        //
        // degree 3:
        // l_last(X) * (z'(X)^2 - z'(X)) = 0
        //
        // which will allow the last value to be zero to
        // ensure the argument is perfectly complete.

        // There are constraints of degree 3 regardless of the
        // number of columns involved.
        3
    }

    pub(crate) fn add_column(&mut self, column: Column<Any>) {
        if !self.columns.contains(&column) {
            self.columns.push(column);
        }
    }

    #[allow(dead_code)]
    pub(crate) fn get_columns(&self) -> Vec<Column<Any>> {
        self.columns.clone()
    }

    pub(crate) fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        assert!(self.columns.len() <= u32::max_value() as usize);
        writer.write_all(&(self.columns.len() as u32).to_le_bytes())?;
        for c in &self.columns {
            writer.write_all(&(c.index() as u64).to_le_bytes())?;
            c.column_type().write(writer)?;
        }
        Ok(())
    }

    pub(crate) fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        let count = {
            let mut data = [0u8; 4];
            reader.read_exact(&mut data[..])?;
            u32::from_le_bytes(data)
        };
        let columns: Vec<_> = (0..count)
            .map(|_| {
                let index = {
                    let mut data = [0u8; 8];
                    reader.read_exact(&mut data[..])?;
                    u64::from_le_bytes(data) as usize
                };
                let column_type = Any::read(reader)?;
                Ok(Column::new(index, column_type))
            })
            .collect::<Result<_, io::Error>>()?;
        Ok(Self { columns })
    }
}

/// The verifying key for a single permutation argument.
#[derive(Clone, Debug)]
pub(crate) struct VerifyingKey<C: CurveAffine> {
    commitments: Vec<C>,
}
impl<C: CurveAffine> VerifyingKey<C> {
    pub(crate) fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
        assert!(self.commitments.len() <= u32::max_value() as usize);
        writer.write_all(&(self.commitments.len() as u32).to_le_bytes())?;
        for c in &self.commitments {
            writer.write_all(c.to_bytes().as_ref())?;
        }
        Ok(())
    }

    pub(crate) fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        let count = {
            let mut data = [0u8; 4];
            reader.read_exact(&mut data[..])?;
            u32::from_le_bytes(data)
        };
        let commitments: Vec<_> = (0..count)
            .map(|_| C::read(reader))
            .collect::<Result<_, _>>()?;
        Ok(Self { commitments })
    }
}
