# -*- coding: utf-8 -*-
from sage.all import ZZ, log


class Cost:
    """
    Algorithms costs.
    """

    # TODO review this list
    _do_repeat = {
        u"rop": True,
        u"red": True,
        u"babai": True,
        u"babai_op": True,
        u"epsilon": False,
        u"mem": False,
        u"delta": False,
        u"beta": False,
        u"k": False,
        u"D_reg": False,
        u"t": False,
        u"m": True,
        u"d": False,
        u"|v|": False,
        u"amplify": False,
        u"repeat": False,  # we deal with it below
        u"c": False,
        u"b": False,
        u"t1": False,
        u"t2": False,
        u"l": False,
        u"ncod": False,
        u"ntop": False,
        u"ntest": False,
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

    def str(self, keyword_width=None, newline=None, round_bound=2048, compact=False, unicode=True):
        """

        :param keyword_width:  keys are printed with this width
        :param newline:        insert a newline
        :param round_bound:    values beyond this bound are represented as powers of two
        :param compact:        do not add extra whitespace to align entries
        :param unicode:        use unicode to shorten representation

        EXAMPLE::

            sage: from estimator import Cost
            sage: s = Cost({"delta_0":5, "bar":2})
            sage: print(s)
            delta_0: 5, bar: 2

            sage: s = Cost([(u"delta_0", 5), ("bar",2)])
            sage: print(s)
            delta_0: 5, bar: 2

        """
        if unicode:
            unicode_replacements = {"delta": u"δ", "beta": u"β", "epsilon": u"ε"}
        else:
            unicode_replacements = {}

        format_strings = {
            u"beta": u"%s: %4d",
            u"d": u"%s: %4d",
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

            sage: from estimator import Cost
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
        u"""
        Return a report with all costs multiplied by `times`.

        :param d:      a cost estimate
        :param times:  the number of times it should be run
        :param select: toggle which fields ought to be repeated and which shouldn't
        :param lll:    if set amplify lattice reduction times assuming the LLL algorithm suffices and costs ``lll``
        :returns:      a new cost estimate

        We maintain a local dictionary which decides if an entry is multiplied by `times` or not.
        For example, δ would not be multiplied but "#bop" would be. This check is strict such that
        unknown entries raise an error. This is to enforce a decision on whether an entry should be
        multiplied by `times` if the function `report` reports on is called `times` often.

        EXAMPLE::

            sage: from estimator import Param, dual
            sage: n, alpha, q = Param.Regev(128)

            sage: dual(n, alpha, q).repeat(2^10)
                rop:   2^91.1
                  m:   2^18.6
                red:   2^91.1
            delta_0: 1.008631
               beta:      115
                  d:      380
                |v|:  688.951
             repeat:   2^27.0
            epsilon: 0.007812

            sage: dual(n, alpha, q).repeat(1)
                rop:   2^81.1
                  m:      380
                red:   2^81.1
            delta_0: 1.008631
               beta:      115
                  d:      380
                |v|:  688.951
             repeat:   2^17.0
            epsilon: 0.007812

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
        return self.str(unicode=False, compact=True)

    def __repr__(self):
        return self.str(unicode=False, newline=True, keyword_width=12)

    def __unicode__(self):
        return self.str(unicode=True)
