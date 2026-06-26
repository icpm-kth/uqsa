# More examples

Here are some more examples:

- [Sample AKAR4 with Parallel
  Tempering](https://icpm-kth.github.io/uqsa/articles/sampleAKAR4mpi.md)
- [Simulate AKAP79 deterministically and add
  noise](https://icpm-kth.github.io/uqsa/articles/simAKAP79.md)
- [Simulate AKAP79
  stochastically](https://icpm-kth.github.io/uqsa/articles/simAKAP79stochastic.md)

It is also possible to use a likelihood-free method (specifically, ABC)
for deterministic models, if measurement noise is taken into account.
So, if we turn on the `noise` option of the simulator, we can just use
ABC on the deterministic (ODE) reaction network models, obtaining
similar results as with the deterministic examples (where
likelihood-based MCMC methods were used):

- [Sample deterministic AKAR4 with
  ABC](https://icpm-kth.github.io/uqsa/articles/AKAR4.md)
- [Sample deterministic AKAP79 with
  ABC](https://icpm-kth.github.io/uqsa/articles/AKAP79.md)
