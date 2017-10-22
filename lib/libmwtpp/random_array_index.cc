#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
/* This little procedure was derived from the die.cpp example in
the boost random library documentation.   It uses the mt19937
random number generator initialized by the current time to create
a unique random number sequence on each invocation.
Following the advice of the documentation the random number generator
is at file scope.  That way the generator is initialized only once
on startup.

This differs from die in two ways.  A die has 6 outcomes.  Here we make
the range an argument = range argument.   Second, the return is intended to
be used as a C array range meaning the return values are random integers
between 0 and range-1.*/
boost::random::mt19937 gen(std::time(0));
int random_array_index(int range) {

    boost::random::uniform_int_distribution<> dist(1, range);
    /*<< A distribution is a function object.  We generate a random
        number by calling `dist` with the generator. Copied from boost
        documentation.
    >>*/
    int dfort=dist(gen);
    return(dfort-1);
}
