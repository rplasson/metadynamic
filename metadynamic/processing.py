#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2019 by Raphaël Plasson
#
# This file is part of metadynamic
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.

from numpy import array, convolve, ones, sqrt


class Result:
    def __init__(self, data):
        res = array(data)
        self.data = {
            "table": res[:, 0],
            "lendist": res[:, 1],
            "pooldist": res[:, 2],
            "end": res[:, 3],
        }
        datasum = {
            "table": self.frame_sum(self.data["table"]),
            "lendist": self.frame_sum(self.data["lendist"]),
            "pooldist": self.frame_sum(self.data["pooldist"]),
        }
        data_n = {
            "table": datasum["table"].loc["#"],
            "lendist": datasum["lendist"].loc["#"],
            "pooldist": datasum["pooldist"].loc["#"],
        }
        self.datamean = {
            "table": datasum["table"] / data_n["table"],
            "lendist": datasum["lendist"] / data_n["lendist"],
            "pooldist": datasum["pooldist"] / data_n["pooldist"],
        }
        self.datastd = {
            "table": sqrt(self.frame_sum([(data-self.datamean["table"])**2 for data in self.data["table"]])/data_n["table"]),
            "lendist": sqrt(self.frame_sum([(data-self.datamean["lendist"])**2 for data in self.data["lendist"]])/data_n["lendist"]),
            "pooldist": sqrt(self.frame_sum([(data-self.datamean["pooldist"])**2 for data in self.data["pooldist"]])/data_n["pooldist"]),
        }

    def _format(self, name, field, num, mean=None):
        if num is None or num == "m":
            res = self.datamean[name]
        elif num == "s":
            res = self.datastd[name]
        elif num == "+":
            res = self.datamean[name] + self.datastd[name]
        elif num == "-":
            res = self.datamean[name] - self.datastd[name]
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

    @staticmethod
    def frame_sum(datas):
        for i, data in enumerate(datas):
            if i == 0:
                res = data
            else:
                res = res.add(data, fill_value=0)
        return res
