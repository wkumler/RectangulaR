# RectangulaR
A repository dedicated to transforming mass-spectrometry data into rectangular format and applying the new format to solving old problems

### Background

Mass spectrometry data consists of 3 dimensions: mass (m/z), retention time (rt), and intensity. Retention time is recorded regularly, and can be referred to by either time (in seconds) or scan number (as an integer). Intriguingly, mass can also be treated in this same way because of spectrometer's tendency to collect m/z in bins, rather than reporting the most accurate mass possible. This binning produces "swaths" of m/z data that are recorded. However, due to the introduction of noise between swath measurements, it's more difficult to assign a single m/z value to an entire swath.

However, an m/z value is assigned to each swath, the data becomes rectangular rather than ragged. I think this would be helpful - we'll see!

There are a few well-defined steps in this process.

1. Convert .raw files to .mzML and separate positive and negative scan modes

2. Identify swaths in (each?) file and replace the noisy m/z values with the median value

3. Perform datum-by-datum blank subtraction

4. Carry on with normal processing (centroiding, peak picking, etc.)
