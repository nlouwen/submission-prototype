""" Collection of custom errors """


class ReferenceNotFound(Exception):
    """Error to signify a reference was not found"""

    def __init__(self, reference) -> None:
        self.message = f"Citation not found: {reference}"
        super().__init__(self.message)
