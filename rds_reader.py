"""Minimal reader for R serialization (.rds / RData, XDR format, versions 2 and 3).

Supports the object types that appear in saved data frames / plain lists:
vectors (logical, integer, double, complex, string, raw), pairlists, generic
vectors (lists), symbols, attributes, NULL, S4 shells, environments (as empty
dicts), ALTREP compact sequences / wrappers / deferred strings, and REFSXP.
Closures / bytecode / promises are read structurally but returned as opaque
placeholders.  Data frames are converted to pandas DataFrames, factors to
pandas Categorical, logicals to nullable booleans.
"""
import gzip, bz2, lzma, struct
import numpy as np
import pandas as pd

NA_INT = -2147483648


class RObj:
    """Generic R object: value + attributes (dict)."""
    __slots__ = ("kind", "value", "attrs", "tag")

    def __init__(self, kind, value=None, attrs=None, tag=None):
        self.kind, self.value, self.attrs, self.tag = kind, value, attrs or {}, tag

    def __repr__(self):
        return f"RObj({self.kind}, attrs={list(self.attrs)})"


class Reader:
    def __init__(self, data: bytes):
        self.buf = data
        self.pos = 0
        self.refs = []

    # ---- primitives ----
    def _int(self):
        v = struct.unpack(">i", self.buf[self.pos:self.pos + 4])[0]
        self.pos += 4
        return v

    def _ints(self, n):
        a = np.frombuffer(self.buf, dtype=">i4", count=n, offset=self.pos).astype(np.int64)
        self.pos += 4 * n
        return a

    def _doubles(self, n):
        a = np.frombuffer(self.buf, dtype=">f8", count=n, offset=self.pos).astype(np.float64)
        self.pos += 8 * n
        return a

    def _bytes(self, n):
        b = self.buf[self.pos:self.pos + n]
        self.pos += n
        return b

    # ---- header ----
    def read(self):
        fmt = self._bytes(2)
        if fmt != b"X\n":
            raise ValueError("only XDR format supported, got %r" % fmt)
        version = self._int()
        self._int()  # writer R version
        self._int()  # min reader version
        if version == 3:
            nelen = self._int()
            self._bytes(nelen)  # native encoding
        return self.item()

    # ---- items ----
    def item(self):
        flags = self._int()
        typ = flags & 0xFF
        levels = flags >> 12
        is_obj = (flags >> 8) & 1
        has_attr = (flags >> 9) & 1
        has_tag = (flags >> 10) & 1

        if typ == 254:  # NILVALUE
            return None
        if typ == 255:  # REFSXP
            idx = flags >> 8
            if idx == 0:
                idx = self._int()
            return self.refs[idx - 1]
        if typ in (242, 241, 243, 253, 252, 251):
            return RObj({242: "emptyenv", 241: "baseenv", 243: "globalenv",
                         253: "unbound", 252: "missingarg", 251: "basenamespace"}[typ])
        if typ in (247, 246, 245):  # PERSISTSXP, PACKAGESXP, NAMESPACESXP
            self._int()  # 0
            n = self._int()
            names = [self.item() for _ in range(n)]
            o = RObj({247: "persist", 246: "package", 245: "namespace"}[typ], names)
            self.refs.append(o)
            return o
        if typ == 238:  # ALTREP
            info = self.item()
            state = self.item()
            attr = self.item()
            return self._altrep(info, state, attr)
        if typ in (240, 239):  # ATTRLISTSXP / ATTRLANGSXP
            typ = 2 if typ == 240 else 6
            has_attr = 1

        if typ == 1:  # SYMSXP
            name = self.item()
            o = RObj("symbol", name.value if isinstance(name, RObj) else name)
            self.refs.append(o)
            return o

        if typ in (2, 6, 17, 3, 5):  # pairlist-like
            attrs = self.item() if has_attr else None
            tag = self.item() if has_tag else None
            car = self.item()
            cdr = self.item()
            kind = {2: "pairlist", 6: "lang", 17: "dots", 3: "closure", 5: "promise"}[typ]
            o = RObj(kind, (car, cdr), tag=tag)
            if attrs is not None:
                o.attrs = pairlist_to_dict(attrs)
            return o

        if typ == 4:  # ENVSXP
            self._int()  # locked
            o = RObj("env", {})
            self.refs.append(o)
            enclos = self.item(); frame = self.item(); hashtab = self.item(); attrib = self.item()
            return o

        if typ == 9:  # CHARSXP
            n = self._int()
            if n == -1:
                return None  # NA_character_
            b = self._bytes(n)
            enc = "latin-1" if (levels & 4) else "utf-8"
            try:
                return b.decode(enc)
            except UnicodeDecodeError:
                return b.decode("latin-1")

        if typ in (10, 13):  # LGLSXP, INTSXP
            n = self._int()
            v = self._ints(n)
            o = RObj("logical" if typ == 10 else "integer", v)
        elif typ == 14:
            n = self._int()
            o = RObj("double", self._doubles(n))
        elif typ == 15:
            n = self._int()
            d = self._doubles(2 * n)
            o = RObj("complex", d[0::2] + 1j * d[1::2])
        elif typ == 16:
            n = self._int()
            o = RObj("string", [self.item() for _ in range(n)])
        elif typ in (19, 20):
            n = self._int()
            o = RObj("list", [self.item() for _ in range(n)])
        elif typ == 24:
            n = self._int()
            o = RObj("raw", self._bytes(n))
        elif typ == 25:
            o = RObj("s4", None)
        elif typ in (7, 8):
            n = self._int()
            o = RObj("builtin", self._bytes(n).decode("latin-1"))
        elif typ == 21:  # BCODESXP -- read structurally
            o = self._bytecode()
        elif typ in (11, 12):  # EXTPTR / WEAKREF
            o = RObj("extptr", None)
            self.refs.append(o)
            if typ == 11:
                self.item(); self.item()
        else:
            raise ValueError(f"unsupported SEXPTYPE {typ} at {self.pos}")

        if has_attr:
            o.attrs = pairlist_to_dict(self.item())
        return o

    def _bytecode(self):
        nreps = self._int()
        reps = [None] * nreps
        return self._bc1(reps)

    def _bc1(self, reps):
        code = self.item()  # INTSXP
        n = self._int()
        consts = []
        for _ in range(n):
            t = self._int()
            if t == 21:
                consts.append(self._bc1(reps))
            elif t in (240, 239, 2, 6):
                consts.append(self._bclang(t, reps))
            else:
                consts.append(self.item())
        return RObj("bytecode", None)

    def _bclang(self, t, reps):
        if t in (241, 240, 239):  # BCREPREF / BCREPDEF / ...
            pass
        if t == 241:  # BCREPREF
            self._int(); return None
        if t == 240:  # BCREPDEF
            self._int(); t = self._int()
        if t in (6, 2, 239, 240):
            if t in (239, 240):
                self.item()  # attr
            self.item()  # tag
            self._bclang(self._int(), reps)  # car
            self._bclang(self._int(), reps)  # cdr
            return None
        return self.item()

    def _altrep(self, info, state, attr):
        cls = info.value[0].value if isinstance(info.value[0], RObj) else str(info.value[0])
        if cls in ("compact_intseq", "compact_realseq"):
            n, start, step = state.value
            arr = start + step * np.arange(int(n))
            o = RObj("integer" if cls == "compact_intseq" else "double",
                     arr.astype(np.int64) if cls == "compact_intseq" else arr.astype(float))
        elif cls == "deferred_string":
            src = state.value[0] if isinstance(state, RObj) and state.kind == "pairlist" else state
            vals = src.value
            o = RObj("string", [None if (src.kind == "integer" and v == NA_INT) or (src.kind == "double" and np.isnan(v))
                                else (str(int(v)) if src.kind == "integer" else repr(float(v))) for v in vals])
        elif cls.startswith("wrap_"):
            o = state.value[0]  # pairlist: car is the wrapped object
        else:
            raise ValueError(f"unsupported ALTREP class {cls}")
        if attr is not None:
            o.attrs = {**o.attrs, **pairlist_to_dict(attr)}
        return o


def pairlist_to_dict(pl):
    d = {}
    while isinstance(pl, RObj) and pl.kind in ("pairlist", "lang"):
        car, cdr = pl.value
        key = pl.tag.value if isinstance(pl.tag, RObj) else pl.tag
        d[key] = car
        pl = cdr
    return d


# ---- conversion to Python/pandas ----
def to_py(o):
    """Convert an RObj tree into Python/pandas structures."""
    if o is None:
        return None
    if not isinstance(o, RObj):
        return o
    cls = o.attrs.get("class")
    cls = cls.value if isinstance(cls, RObj) else []
    if o.kind == "list":
        names = o.attrs.get("names")
        names = names.value if isinstance(names, RObj) else None
        if "data.frame" in cls:
            cols = {nm: to_py(col) for nm, col in zip(names, o.value)}
            df = pd.DataFrame(cols)
            return df
        if names is not None:
            return {nm: to_py(v) for nm, v in zip(names, o.value)}
        return [to_py(v) for v in o.value]
    if o.kind == "integer":
        if "factor" in cls:
            levels = o.attrs["levels"].value
            codes = np.where(o.value == NA_INT, -1, o.value - 1)
            return pd.Categorical.from_codes(codes, categories=levels)
        v = o.value
        if (v == NA_INT).any():
            return pd.array(np.where(v == NA_INT, 0, v), dtype="Int64").where(v != NA_INT, pd.NA) if False else \
                pd.Series(v, dtype="int64").mask(v == NA_INT).astype("Int64").to_numpy()
        return v
    if o.kind == "logical":
        v = o.value
        out = pd.array([None if x == NA_INT else bool(x) for x in v], dtype="boolean")
        return out
    if o.kind == "double":
        return o.value
    if o.kind == "string":
        return pd.array(o.value, dtype="string") if any(x is None for x in o.value) else np.array(o.value, dtype=object)
    if o.kind == "pairlist":
        return {k: to_py(v) for k, v in pairlist_to_dict(o).items()}
    if o.kind == "symbol":
        return o.value
    return o


def _decompress(raw: bytes) -> bytes:
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw)
    if raw[:3] == b"BZh":
        return bz2.decompress(raw)
    if raw[:6] == b"\xfd7zXZ\x00":
        return lzma.decompress(raw)
    return raw


def read_rds(path):
    with open(path, "rb") as f:
        raw = _decompress(f.read())
    return to_py(Reader(raw).read())


def read_rdata(path):
    """Read an RData/.rda file (saved with save()); returns dict name->object."""
    with open(path, "rb") as f:
        raw = _decompress(f.read())
    if raw[:5] != b"RDX2\n" and raw[:5] != b"RDX3\n":
        raise ValueError("not an RData file")
    r = Reader(raw[5:])
    obj = r.read()
    out = {}
    pl = obj
    while isinstance(pl, RObj) and pl.kind == "pairlist":
        car, cdr = pl.value
        out[pl.tag.value] = to_py(car)
        pl = cdr
    return out


if __name__ == "__main__":
    import sys
    for p in sys.argv[1:]:
        o = read_rdata(p) if p.lower().endswith((".rdata", ".rda")) else read_rds(p)
        if isinstance(o, pd.DataFrame):
            print(p, o.shape); print(o.dtypes.to_string()); print(o.head())
        else:
            print(p, type(o), o if not isinstance(o, dict) else {k: (v.shape if hasattr(v, 'shape') else v) for k, v in o.items()})
