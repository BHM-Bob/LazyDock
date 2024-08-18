import uuid


def uuid4():
    """return a random 4-character string"""
    return uuid.uuid4().hex[:4]