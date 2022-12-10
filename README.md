# Expanding and Fingerprinting a Repeating Earthquake Catalog on the Cascades Volcanoes

We look at the REDPy (<a href="github.com/ahotovec/REDPy"> Repeating Earthquake Detector in Python - Alicia Hotovec-Ellis software</a>) catalog on five Cascades Volcanoes (Mt Rainier, Mt St Helens, Mt Hood, Mt Baker, Newberry Volcano). Because REDPy catalogs only have cluster ID and guessed origin time, we aim to do 3 things:


1. Backfill the REDPy catalog with template matching from 2002-2009
2. Locate the events 
3. Fingerprint events using clustering

We use eqcorrscan (<a href="https://www.dropbox.com/s/rscu5odvn1bbr2s/Chamberlain18.pdf?dl=0"> Chamberlain et al, 2018</a>) and the Fast Matched Filter (<a href="https://doi.org/10.1785/0220170181"> Beauce and Frank, 2018 < /a>) to run template matching on GPUs.

<h2>Understanding the Repo</h2>
<table>
  <tr>
    <th>Filename</th>
    <th>Directory</th>
    <th>Description</th>
  </tr>
  <tr>
    <td>Volcano_Metadata.csv</td>
    <td>csv_catalogs/</td>
    <td>Contains metadata for stations used on all volcanoes, created in make_volcano_metadata</td>
  </tr>
  <tr>
    <td>make_volcano_metadata.ipynb</td>
    <td>notebooks/</td>
    <td>Queries IRIS for station metadata in a certain radius of volcano center and sorts through to find useful stations, makes Volcano_Metadata.csv</td>
  </tr>
</table>