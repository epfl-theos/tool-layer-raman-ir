# Layer-Raman-IR online tool
A tool to compute and visualize Raman-active and Infrared-active modes in layered materials

[![Actions Status](https://github.com/epfl-theos/tool-layer-raman-ir/workflows/Continuous%20integration/badge.svg)](https://github.com/epfl-theos/tool-layer-raman-ir/actions)

## About the tool

This tool allows users to upload the bulk crystal structure of a layered material in a number of common formats
(or to choose from a few examples) and, after detection of bonds and of 2D layers, it determines the symmetry
of the inter-layer force-constant matrices and the corresponding optical-activity fan diagram.

The output page that is displayed includes relevant information on the structure (interactive visualizations
of the bulk multilayer and of each layer, information on the coincidence operation and of the symmetry properties of the multilayer),
and shows the independent components of the force-constant matrices.
For the latter, random values are proposed, chosen so as to mimick their typical ratios.
These can be changed interactively, and the tool then computes in real-time the corresponding fan diagram,
including the optical activity for infrared and Raman spectroscopy
(and, for Raman modes, it also shows if the mode would be detectable in a back-scattering geometry).

Note: the tool only works for structures that satisfy the maximum-degree-of-order hypotheses of order-disorder polytypes,
and in particular that the structure is composed by a stacking of the same 2D layer,
with the same coincidence relationship bringing any layer onto the next one.
If any of the conditions is not satisfied, the tool will display a message informing that the structure does not satisfy the assumptions,
and so the symmetry analysis cannot be applied.

## Online version
This tool is deployed on the Materials Cloud "Tools" section [here](https://materialscloud.org/work/tools/layer-raman-ir), so you can use it without need of installation.

## How to cite
If you use this tool, please cite the following work:

* G. Pizzi, S. Milana, A. C. Ferrari, N. Marzari, M. Gibertini, *Shear and Breathing Modes of Layered Materials*, ACS Nano (2021), [doi:10.1021/acsnano.0c10672](https://doi.org/10.1021/acsnano.0c10672).

You might also want to cite the spglib, ASE and pymatgen libraries that are used internally by the tool.

## How to deploy on your computer
1. Install [Docker](https://www.docker.com)
2. Clone this repository
3. Run `./admin-tools/build-and-run.sh`. This will build the Docker image, start the docker container, and open a browser to the correct URL.
   If the browser does not open automatically, connect to http://localhost:8091

If you want to check the Apache logs, you can run `./admin-tools/get-apache-logs.sh --reload` to see the logs of the tool in real time while you use it.

## Acknowledgements
This tool uses the ASE and pymatgen libraries for structure manipulation, and spglib for symmetry detection.
The tool is based upon the [tools-barebone framework](https://github.com/materialscloud-org/tools-barebone) developed by the Materials Cloud team.

We acknowledge funding from the [MARVEL National Centre of Competence in Research](https://nccr-marvel.ch) of the
Swiss National Science Foundation (SNSF), the
[European Centre of Excellence MaX "Materials design at the Exascale"](http://www.max-centre.eu),
the [swissuniversities P-5 "Materials Cloud" project](https://www.materialscloud.org/swissuniversities),
the [Graphene Flagship](http://graphene-flagship.eu),
the RC Grant Hetero2D, the EPSRC Grants EP/509K01711X/1, EP/K017144/1, EP/N010345/1, EP/M507799/5101 and EP/L016087/1,
the Italian Ministry for University and Research through the Levi-Montalcini program,
and from SNSF through the Ambizione program.
We thank L. Talirz for support during the deployment of this tool on Materials Cloud.
