# OpenSimplexNoise

A C++ port of Open Simplex Noise. The port is based on a implementation in Java by Kurt Spencer.
The implementation includes 2D, 3D and 4D noise.

## Example
Example animation made by the noise.

![](/Examples/animation.gif)

## Usage
Add the .h and .cpp file to your project.
```cpp
OpenSimplexNoise::Noise noise
double value = noise.eval(0.13, 0.31) // For 2D noise
```
