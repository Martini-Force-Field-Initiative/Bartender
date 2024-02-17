
# `write_bartender_inp.py`
A small scripts that converts an AA-to-CG mapping (`ndx`) and CG topology (`itp`) to a Bartender input file (`inp`).

How to use the `write_bartender_inp.py` script:
```
python write_bartender_inp.py --ndx inputs/BENZ_oplsaaTOcg_cgbuilder.ndx         --itp inputs/BENZ_cog.itp --out BENZ_bartender.inp  
python write_bartender_inp.py --ndx inputs/NDMBI_oplsaaTOcg_cgbuilder.ndx        --itp inputs/NDMBI.itp    --out NDMBI_bartender.inp  
python write_bartender_inp.py --ndx inputs/TOLU_oplsaaTOcg_cgbuilder_refined.ndx --itp inputs/TOLU.itp     --out TOLU_bartender.inp  
```

