# Concentration to Count Conversion

This function returns the factor f for converting a concentration to
particle numbers (counts) based on the units specified in the SBtab file
and the volume of the system encoded as LV (Avogadro's constant ×
Volume).

## Usage

``` r
ccc(sb)
```

## Arguments

- sb:

  SBtab list of data.frames with one of the members called 'Compound',
  with a mandatory "!Unit" column.

- LV:

  Avogadro's constant multiplied by system volume

## Value

f a conversion factor: n \<- f\*x, where n is a particle count and x a
concentration in the specified volume.

## Details

The implicit assumption is that the LV constant the user passes to this
function has a compatible unit of volume to the SBtab unit of
concentration. LV should be in litres and the concentrations also in
litres up to SI prefixes. All of these are OK: M, nmol/l, mol/dl,
kilomol/megalitre.

If the concentrations are 'something per cubic metre', then LV has to be
in m³ as well. The unit of LV is not checked in any way. But the unit of
concentration will be parsed an converted to mol/l.
