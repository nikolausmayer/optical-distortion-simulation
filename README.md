
# optical-distortion-simulation
Approximate simulation and inversion of optical camera distortion, using binary search. This is old crappy code (it only covers the radial part of a radial-tangential distortion model), but it works ✌!

# About 
The application will take your image, distort it, and then *un*distort the intermediate result to arrive back at the beginning -- except that a whole lot of information is lost along the way.

The final image is an approximation of what an undistorted image of a real camera looks like, given the specified distortion parameters. This code was used to generate training data for [a stereo-vision deep neural network](https://lmb.informatik.uni-freiburg.de/Publications/2019/MB19/) in the ["TrimBot" EU Horizon 2020 research project](https://lmb.informatik.uni-freiburg.de/research/funded_projects/eu_trimbot/), also also for [another paper](https://lmb.informatik.uni-freiburg.de/Publications/2018/MIFDB18/) and finally in [my PhD thesis](https://freidok.uni-freiburg.de/data/166944).

# Usage
Put [CImg.h](http://cimg.eu/) into `src`, then `make` and run the resulting `generic_executable`. 

