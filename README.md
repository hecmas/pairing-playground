# Pairing Playground

This repo contains various implementations of the evolution of pairings as explained in the book [Pairings for Beginners](https://static1.squarespace.com/static/5fdbb09f31d71c1227082339/t/5ff394720493bd28278889c6/1609798774687/PairingsForBeginners.pdf) by Craig Costello.

In particular:
- `weil_tate_pairing_naive.sage` implements Miller's algorithm for the Weil and the Tate pairings "as prescribed" (algorithm 5.1 in the book). That is, there are no optimizations over the initial definitions of both Weil and Tate pairings beyond the usage of Miller's double-and-add algorithm. The file `divisors.sage` contains the `Divisors` class useful for playing around with divisors.
- `tate_pairing.sage` contains the BKLS-GHS version of Miller's algorithm for the Tate pairing (algorithm 7.1 in the book). In this version, divisors are not needed anymore and the irrelevant factors technique has been used.
- `ate_pairing.sage` contains the BKLS-GHS version of Miller's algorithm for the ate pairing (algorithm 7.2 in the book). The ate pairing is the Tate pairing with a much shorter Miller's loop and domains reversed.

The file `tools.sage` contains various methods shared among all the previous files.

## TODOs

- [x] Divisors.
- [x] Weil pairing.
- [x] Tate pairing.
- [x] Ate pairing.
- [ ] Optimal pairings. In particular, the optimal ate pairing.
- [ ] Final exponentiation optimizations.
- [ ] Twists.
- [ ] Endomorphisms.
- [ ] Towered extension fields.
- [ ] Projective coordinates.