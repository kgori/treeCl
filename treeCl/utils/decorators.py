def lazyprop(fn):
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazyprop(self):
        if hasattr(self, attr_name):
            return getattr(self, attr_name)
        else:
            value = fn(self)
            if value:
                setattr(self, attr_name, value)
                return value

    return _lazyprop
