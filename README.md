# Layer-Raman online tool
A tool to compute and visualize Raman-active and Infrared-active modes in layered materials

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
This tool is deployed on the Materials Cloud "Tools" section [here](https://layer-raman.materialscloud.io), so you can use it without need of installation.

## How to cite
If you use this tool, please cite the following work:

* G. Pizzi, S. Milana, A. C. Ferrari, N. Marzari, M. Gibertini, *Shear and breathing modes of all layered materials*, to be submitted (2020).

You might also want to cite the spglib, ASE and pymatgen libraries that are used internally by the tool.

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
