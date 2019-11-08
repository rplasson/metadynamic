from numpy import array, convolve, ones


class Result:
    def __init__(self, data):
        res = array(data)
        self.data = {
            "table": res[:, 0],
            "lendist": res[:, 1],
            "pooldist": res[:, 2],
            "end": res[:, 3],
        }

    def _format(self, name, field, num, mean=None):
        if num is None or num == "m":
            res = self.data[name].mean()
        elif num == "s":
            res = self.data[name].std()
        elif num == "+":
            res = self.data[name].mean() + self.data[name].std()
        elif num == "-":
            res = self.data[name].mean() - self.data[name].std()
        else:
            res = self.data[name][num]
        if field is not None:
            res = res.loc[field]
        if mean:
            res = self.running_mean(res, mean)
        return res

    def table(self, field=None, num=None, mean=None):
        return self._format("table", field, num, mean)

    def lendist(self, field=None, num=None, mean=None):
        return self._format("lendist", field, num, mean)

    def pooldist(self, field=None, num=None, mean=None):
        return self._format("pooldist", field, num, mean)

    def end(self, num=None):
        if num is None:
            return self.data["end"]
        return self.data["end"][num]

    @staticmethod
    def running_mean(data, length):
        """Running mean from
           https://stackoverflow.com/questions/13728392/moving-average-or-running-mean"""
        return convolve(data, ones((length,)) / length, mode="valid")
