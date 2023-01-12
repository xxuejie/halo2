use std::io;

use pasta_curves::{arithmetic::CurveAffine, Fp, Fq};

pub(crate) trait CurveRead: CurveAffine {
    /// Reads a compressed element from the buffer and attempts to parse it
    /// using `from_bytes`.
    fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
        let mut compressed = Self::Repr::default();
        reader.read_exact(compressed.as_mut())?;
        Option::from(Self::from_bytes(&compressed))
            .ok_or_else(|| io::Error::new(io::ErrorKind::Other, "invalid point encoding in proof"))
    }
}

impl<C: CurveAffine> CurveRead for C {}

/// This serves as a temporary solution till `from_raw` can become a
/// generic trait method of scalar elements
pub trait FieldEncoding: Sized {
    fn read<R: io::Read>(reader: &mut R) -> io::Result<Self>;
    fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()>;
}

macro_rules! impl_field_encoding {
    ($name:ident) => {
        impl FieldEncoding for $name {
            fn read<R: io::Read>(reader: &mut R) -> io::Result<Self> {
                let mut data = [0u8; 32];
                reader.read_exact(&mut data[..])?;

                let mut raw = [0u64; 4];
                for i in 0..32 {
                    raw[i / 8] |= (data[i] as u64) << ((i % 8) * 8);
                }
                Ok(Self::from_raw(raw))
            }

            fn write<W: io::Write>(&self, writer: &mut W) -> io::Result<()> {
                let data: [u8; 32] = self.into();
                writer.write_all(&data)
            }
        }
    };
}

impl_field_encoding!(Fp);
impl_field_encoding!(Fq);
