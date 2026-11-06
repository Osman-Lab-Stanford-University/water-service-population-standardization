# Saved-output verification of standardized CDFs

The 41-threshold CDF table was compared with the E23 JSON result register at the CSV's reported precision. E15 estimates the standardized CDFs; E23 consolidates those results.

| Check | Result |
|---|---|
| Five source-specific CDFs × 41 thresholds | 205 of 205 values agreed |
| Piped-minus-groundwater differences | 41 of 41 values agreed |
| Total checked fields | 246, with no differences exceeding the rounding tolerance |
| Maximum absolute difference | 0.00049971475101 percentage points |
| Rounding tolerance | 0.000500001 percentage points |
| Archived E23 JSON register SHA-256 | `27f962002af43cb95643fe95920fc2af14a973abbc985e305ba8c9b87aab9cf6` |
| E23 register manifest | Identity confirmed |
| Separate downstream 25H comparison | 246 matching fields |
| Archived 25H source SHA-256 | `978d0e967d93a0b41df1c2d7b51692305959fa1308e9d4b38b462436881fd5d6` |
| Released CDF CSV SHA-256 | `e2d6c96b0ebf7dd83ca4a49f688fd65838590d45b3b8106c2a53d073ce247f16` |

The tolerance reflects rounding of the CDF values recorded in the aggregate CSV.
