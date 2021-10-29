# -*- coding: utf-8 -*-
from sage.all import ZZ, log


class Cost:
    """
    Algorithms costs.
    """

    # TODO review this list
    _do_repeat = {
        "rop": True,
        "red": True,
        "delta": False,
        "beta": False,
        "eta": False,
        "epsilon": False,
        "m": True,
        "d": False,
        "amplify": False,
        "repeat": False,  # we deal with it below
    }

    def __init__(self, data=None, **kwds):
        """

        :param data: we call ``dict(data)``

        """
        if data is None:
            self.data = dict()
        else:
            self.data = dict(data)

        for k, v in kwds.items():
            self.data[k] = v

    def str(  # noqa C901
        self, keyword_width=None, newline=None, round_bound=2048, compact=False, unicode=True
    ):
        """

        :param keyword_width:  keys are printed with this width
        :param newline:        insert a newline
        :param round_bound:    values beyond this bound are represented as powers of two
        :param compact:        do not add extra whitespace to align entries
        :param unicode:        use unicode to shorten representation

        EXAMPLE::

            sage: from estimator.cost import Cost
            sage: s = Cost({"delta":5, "bar":2})
            sage: print(s)
            δ: 5, bar: 2

            sage: s = Cost([(u"delta", 5), ("bar",2)])
            sage: print(s)
            δ: 5, bar: 2

        """
        if unicode:
            unicode_replacements = {"delta": "δ", "beta": "β", "eta": "η", "epsilon": "ε"}
        else:
            unicode_replacements = {}

        format_strings = {
            "beta": "%s: %4d",
            "d": "%s: %4d",
            "b": "%s: %3d",
            "t1": "%s: %3d",
            "t2": "%s: %3d",
            "l": "%s: %3d",
            "ncod": "%s: %3d",
            "ntop": "%s: %3d",
            "ntest": "%s: %3d",
        }

        d = self.data
        s = []
        for k in d:
            if k == "tag":
                fmt = u"%%%ds" % keyword_width
                if compact:
                    s.append("%s: %s" % (fmt % k, d[k]))
                else:
                    s.append("%s: %8s" % (fmt % k, d[k]))
                continue
            v = d[k]
            kk = unicode_replacements.get(k, k)
            if keyword_width:
                fmt = u"%%%ds" % keyword_width
                kk = fmt % kk
            if not newline and k in format_strings:
                s.append(format_strings[k] % (kk, v))
            elif (
                ZZ(1) / round_bound < v < round_bound
                or v == 0
                or ZZ(-1) / round_bound > v > -round_bound
            ):
                try:
                    if compact:
                        s.append(u"%s: %d" % (kk, ZZ(v)))
                    else:
                        s.append(u"%s: %8d" % (kk, ZZ(v)))
                except TypeError:
                    if v < 2.0 and v >= 0.0:
                        if compact:
                            s.append(u"%s: %.6f" % (kk, v))
                        else:
                            s.append(u"%s: %8.6f" % (kk, v))
                    else:
                        if compact:
                            s.append(u"%s: %.3f" % (kk, v))
                        else:
                            s.append(u"%s: %8.3f" % (kk, v))
            else:
                t = u"%s" % (u"≈" if unicode else "") + u"%s2^%.1f" % (
                    "-" if v < 0 else "",
                    log(abs(v), 2).n(),
                )
                if compact:
                    s.append(u"%s: %s" % (kk, t))
                else:
                    s.append(u"%s: %8s" % (kk, t))
        if not newline:
            if compact:
                return u", ".join(s)
            else:
                return u",  ".join(s)
        else:
            return u"\n".join(s)

    def reorder(self, first):
        """
        Return a new ordered dict from the key:value pairs in dictinonary but reordered such that the
        ``first`` keys come first.

        :param dictionary: input dictionary
        :param first: keys which should come first (in order)

        EXAMPLE::

            sage: from estimator.cost import Cost
            sage: d = Cost([("a",1),("b",2),("c",3)]); d
            a:        1
            b:        2
            c:        3

            sage: d.reorder( ["b","c","a"])
            b:        2
            c:        3
            a:        1
        """
        keys = list(self.data)
        for key in first:
            keys.pop(keys.index(key))
        keys = list(first) + keys
        r = dict()
        for key in keys:
            r[key] = self.data[key]
        return Cost(r)

    def filter(self, keys):
        """
        Return new ordered dictinonary from dictionary restricted to the keys.

        :param dictionary: input dictionary
        :param keys: keys which should be copied (ordered)
        """
        r = dict()
        for key in keys:
            r[key] = self.data[key]
        return Cost(r)

    def repeat(self, times, select=None, lll=None):
        """
        Return a report with all costs multiplied by ``times``.

        :param times:  the number of times it should be run
        :param select: toggle which fields ought to be repeated and which should not
        :param lll:    if set amplify lattice reduction times assuming the LLL algorithm suffices and costs ``lll``
        :returns:      a new cost estimate

        We maintain a local dictionary which decides if an entry is multiplied by ``times`` or not.
        For example, ``δ`` would not be multiplied but ``rop`` would be. This check is strict such that
        unknown entries raise an error. This is to enforce a decision on whether an entry should be
        multiplied by ``times`` if the function ``report`` reports on is called `times` often.

        """

        if lll and self["red"] != self["rop"]:
            raise ValueError("Amplification via LLL was requested but 'red' != 'rop'")

        if select is not None:
            for key in select:
                self._do_repeat[key] = select[key]

        ret = dict()
        for key in self.data:
            try:
                if self._do_repeat[key]:
                    if lll and key in ("red", "rop"):
                        ret[key] = self[key] + times * lll
                    else:
                        ret[key] = times * self[key]
                else:
                    ret[key] = self.data[key]
            except KeyError:
                raise NotImplementedError(
                    u"You found a bug, this function does not know about '%s' but should." % key
                )
        ret[u"repeat"] = times * ret.get("repeat", 1)
        return Cost(ret)

    def __rmul__(self, times):
        return self.repeat(times)

    def combine(self, right, base=None):
        """Combine ``left`` and ``right``.

        :param left: cost dictionary
        :param right: cost dictionary
        :param base: add entries to ``base``

        """
        if base is None:
            cost = Cost()
        else:
            cost = base
        for key in self.data:
            cost[key] = self.data[key]
        for key in right:
            cost[key] = right.data[key]
        return Cost(cost)

    def __add__(self, other):
        return self.combine(self, other)

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value

    def __iter__(self):
        return iter(self.data)

    def values(self):
        return self.data.values()

    def __str__(self):
        return self.str(unicode=True, compact=True)

    def __repr__(self):
        return self.str(unicode=True, newline=True, keyword_width=12)

    def __unicode__(self):
        return self.str(unicode=True)
