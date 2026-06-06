# helpers.py

import math


# =========================================
# SAFE VARIABLE VALUE
# =========================================

def var_value_safe(
    var,
    default=0.0
):

    try:

        if var is None:

            return float(default)

        v = getattr(
            var,
            'X',
            None
        )

        if v is None:

            v = getattr(
                var,
                'x',
                None
            )

        return (

            float(v)

            if v is not None

            else float(default)

        )

    except Exception:

        return float(default)


# =========================================
# SAFE DICTIONARY GET
# =========================================

def safe_get(
    dictionary,
    key,
    default=0.0
):

    try:

        return dictionary.get(
            key,
            default
        )

    except Exception:

        return default


# =========================================
# MINUTE BIN INTEGRATION
# =========================================

def integrate_bin(
    arrivals,
    t1_min,
    t2_min
):

    if t2_min <= t1_min:

        return 0.0

    total = 0.0

    start_min = int(
        math.floor(t1_min)
    )

    end_min = int(
        math.floor(t2_min)
    )

    # =====================================
    # Same Minute
    # =====================================

    if start_min == end_min:

        frac = (
            t2_min
            - t1_min
        )

        if (

            0 <= start_min

            < len(arrivals)

        ):

            return (

                arrivals[start_min]
                * frac

            )

        return 0.0

    # =====================================
    # First Partial Minute
    # =====================================

    first_frac = (

        start_min + 1
        - t1_min

    )

    if (

        0 <= start_min

        < len(arrivals)

    ):

        total += (

            arrivals[start_min]
            * first_frac

        )

    # =====================================
    # Full Minutes
    # =====================================

    for m in range(

        start_min + 1,
        end_min

    ):

        if (

            0 <= m

            < len(arrivals)

        ):

            total += arrivals[m]

    # =====================================
    # Last Partial Minute
    # =====================================

    last_frac = (

        t2_min
        - end_min

    )

    if (

        0 <= end_min

        < len(arrivals)

    ):

        total += (

            arrivals[end_min]
            * last_frac

        )

    return total


# =========================================
# CONVERT SECONDS TO MINUTES
# =========================================

def seconds_to_minutes(
    seconds
):

    try:

        return float(seconds) / 60.0

    except Exception:

        return 0.0


# =========================================
# CONVERT MINUTES TO SECONDS
# =========================================

def minutes_to_seconds(
    minutes
):

    try:

        return float(minutes) * 60.0

    except Exception:

        return 0.0


# =========================================
# CLIP VALUE
# =========================================

def clip(
    value,
    lower,
    upper
):

    return max(
        lower,
        min(
            value,
            upper
        )
    )


# =========================================
# AVERAGE ABSOLUTE DIFFERENCE
# =========================================

def average_absolute_difference(
    dict1,
    dict2
):

    keys = set(
        dict1.keys()
    ).intersection(
        dict2.keys()
    )

    if len(keys) == 0:

        return 0.0

    diff = sum(

        abs(
            dict1[k]
            - dict2[k]
        )

        for k in keys

    )

    return diff / len(keys)


# =========================================
# POSITIVE VALUE FILTER
# =========================================

def positive_entries(
    dictionary
):

    return {

        k: v

        for k, v
        in dictionary.items()

        if v > 0

    }


# =========================================
# SORT VARIABLE NAME
# =========================================

def variable_sort_key(var):

    try:

        name = var.VarName

        inside = name[
            name.find('[') + 1:
            name.rfind(']')
        ]

        parts = [

            p.strip()

            for p in inside.split(',')

        ]

        return tuple(

            float(p)

            if '.' in p

            else int(p)

            for p in parts

        )

    except Exception:

        return (var.VarName,)