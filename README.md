# Caspase Model

This is a modified version of the caspase model based on [Albeck2008](#Albeck2008).
The implemented modifications have been published in [Corbat2018](#Corbat2018).

The model is implemented twice:

- with [PySB](https://pysb.org), in `caspase_model.models`, building on [EARM](https://github.com/acorbat/earm),
- with [SimBio](https://github.com/dyscolab/simbio), in `caspase_model.simbio_model`.

Tests check that both implementations give the same results.

## SimBio implementation

The SimBio models are composed from subsystem modules (`caspase_model.simbio_model.modules`):

| Module | Description |
| --- | --- |
| `ReceptorToBid` | Receptor ligation to Bid truncation |
| `MOMP11b` to `MOMP11e` | Mitochondrial outer membrane permeabilization variants (Albeck2008, Figure 11) |
| `CaspaseCascade` | Smac and CytoC activation, effector caspases and PARP cleavage |
| `Apoptosome` | Apoptosome formation (Albeck2008) |
| `DoubleApoptosome` | Two-step apoptosome activation (ARM) |
| `CytosolicBcl2` | Cytosolic Bcl2 sequestering tBid |
| `IntrinsicStimulus` | Intrinsic stimulus truncating Bid |
| `Biosensor` | Homo-FRET caspase biosensor (Corbat2018) |

Each module owns the species it introduces,
and a model links the shared ones:

```python
from simbio import System, Variable, initial

from caspase_model.simbio_model.modules import (
    MOMP11e,
    Apoptosome,
    Biosensor,
    CaspaseCascade,
    ReceptorToBid,
)


class Model(System):
    L: Variable = initial(default=1e3)

    momp = MOMP11e()
    rec_to_bid = ReceptorToBid(L=L, Bid_U=momp.Bid_U, Bid_T=momp.Bid_T)
    caspases = CaspaseCascade(
        Smac_C=momp.Smac_C,
        CytoC_C=momp.CytoC_C,
        C8_pro=rec_to_bid.C8_pro,
        C8_A=rec_to_bid.C8_A,
    )
    apoptosome = Apoptosome(
        CytoC_A=caspases.CytoC_A,
        C3_pro=caspases.C3_pro,
        C3_A=caspases.C3_A,
        XIAP=caspases.XIAP,
    )
    sCas3 = Biosensor(enzyme=caspases.C3_A, kf=5e-7, kr=1e-2)
```

The published models are in `caspase_model.simbio_model.models`:
`Albeck11b` to `Albeck11f`, `AlbeckAsMatlab`, `Corbat2018` and `ARM`.
They are stimulated extrinsically by default:

```python
from simbio import Simulator

from caspase_model.simbio_model import ARM, intrinsic

extrinsic_result = Simulator(ARM).solve(save_at=range(0, 20_000, 100))
intrinsic_result = Simulator(intrinsic(ARM)).solve(save_at=range(0, 20_000, 100))
```

Observables summing each caspase state over all its complexes,
as `caspase_model.shared.observe_caspases` does for PySB,
are available as simulation outputs:

```python
from caspase_model.simbio_model import observables

result = (
    Simulator(ARM)
    .with_transform(observables(ARM))
    .solve(save_at=range(0, 20_000, 100))
)
```

## References

<a name="Albeck2008">Albeck2008</a>: [Albeck JG, Burke JM, Spencer SL, Lauffenburger DA, Sorger PK (2008) Modeling a Snap-Action, Variable-Delay Switch Controlling Extrinsic Cell Death. PLOS Biology 6(12): e299.](https://doi.org/10.1371/journal.pbio.0060299)

<a name="Corbat2018">Corbat2018</a>: [Corbat AA, Schuermann KC, Liguzinski P , Radon Y, Bastiaens PIH, Verveer PJ, Grecco HE, (2018) Co-imaging extrinsic, intrinsic and effector caspase activity by fluorescence anisotropy microscopy, Redox Biology, Volume 19, Pages 210-217](https://doi.org/10.1016/j.redox.2018.07.023)

## Development

This project is managed by [pixi](https://pixi.sh).
PySB and BioNetGen are installed from the `alubbock` conda channel.

Run the tests with:

```sh
pixi run test
```

Code style is enforced through [prek](https://prek.j178.dev) hooks.
Install them with:

```sh
pixi run prek install
```

or run them on all files with:

```sh
pixi run lint
```
