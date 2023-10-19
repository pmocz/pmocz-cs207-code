Philip Mocz, Xinyi Guo
Harvard University
CS 207, Spring 2014

== Description ==

An undirected graph class.


== HW notes ==

HW1:
  - please run my designed predicate in question #4 on the large data set (it is designed for that one):
      ./subgraph data/large.nodes  data/large.tets 

  - please run the following to generate the `cool image' for question #5:
      ./hw1q5_cool_image data/large.nodes  data/large.tets 
    or take a look at the screenshot
      hw1q5_cool_image.png

HW4:
   run shallow_water.cpp, which uses the our mesh class (Mesh.hpp)




== Files ==

primes.cpp
  calculates number of primes less than input n. Uses memoization to store a vector of divisor candidates, which gets reduced to a list of primes as is_prime() is called for n=2,3,4,... and hence the scaling is better than n^(3/2)

Graph.hpp
  basic graph class
