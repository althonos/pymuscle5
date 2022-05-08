from . import (
    test_muscle,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_muscle))
    return suite
