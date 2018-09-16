"""
Errors raised by multiprot modules
"""

class RanchError(Exception):
    """ Base class for exceptions."""
    pass

class InputError(RanchError):
    """ Exception raised for errors when the input is invalid"""
    pass

class MatchError(RanchError):
    """ Exception raised for errors when comparing two sequences"""
    pass

