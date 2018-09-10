"""
Errors raised by ranch script
"""

class RanchError(Exception):
    """ Base class for exceptions."""
    pass

class MatchError(RanchError):
    """ Exception raised for errors when comparing two sequences"""
    pass


#### Do I need tests here?