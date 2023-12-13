# Expanding and Fingerprinting a Repeating Earthquake Catalog on the Cascades Volcanoes

We look at the REDPy <a href="https://github.com/ahotovec/REDPy">(Repeating Earthquake Detector in Python - Alicia Hotovec-Ellis software)</a> catalog on five Cascades Volcanoes (Mt Rainier, Mt St Helens, Mt Hood, Mt Baker, Newberry Volcano). Because REDPy catalogs only have cluster ID and guessed origin time, we aim to do 3 things:


1. Backfill the REDPy catalog with template matching from 2002-2009
2. Locate the events 
3. Fingerprint events using clustering

We use eqcorrscan <a href="https://www.dropbox.com/s/rscu5odvn1bbr2s/Chamberlain18.pdf?dl=0">(Chamberlain et al, 2018)</a>  and the Fast Matched Filter <a href="https://doi.org/10.1785/0220170181">(Beauce and Frank, 2018)</a> to run template matching on GPUs.

<h2>Methods/Process</h2>
We select stations with make_volcano_metadata.ipynb, generate templates in make_templates.ipynb, locate templates in template_locations.ipynb, detect with templates in find_detections.ipynb, make events in find_events.ipynb, and do post processing in remove_redpy_overlap.ipynb.
<img src="methods.jpg" alt="Methods Flowchart">

<h2>Understanding the Repo</h2>
<table>
  <tr>
    <th>Directory</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>catalogs/</td>
    <td>REDPy catalogs as CSVs with column titles. A subdirectory has the REDPy as they were downloaded (.txt) on June 19, 2022.</td>
  </tr>
  <tr>
    <td>notebooks/</td>
    <td>Jupyter Notebooks of code. Directory is subdivided into backfilling, clustering, location, plots, and tutorials.</td>
  </tr>
  <tr>
    <td>scripts/</td>
    <td>script versions of notebooks, functions (such as mbf_elep_func.py), and config files.</td>
  </tr>
</table>
