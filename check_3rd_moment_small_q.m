
/*
    check_fq_moments_magma.m

    Magma translation of the fast Python DP checker for the theorem:

        For every k-subset I of GF(q), check whether all triples
        (sum x, sum x^2, sum x^3) are attained.

    Default range:
        q <= 73, q >= 17, characteristic >= 5, and only k <= q/2
        by complement symmetry.

    Important output:
        The revised statement with q >= 17 fails for q=17, k=6,7.
        It should pass for the other q <= 73 in characteristic >= 5.

    Run in Magma by:
        magma check_fq_moments_magma.m
*/

SetVerbose("User1", 0);

// -----------------------------------------------------------------------------
// User parameters
// -----------------------------------------------------------------------------

MaxQ := 50;
MinQ := 17;
MinChar := 5;
MinK := 6;

UseSymmetry := true;            // if true, only check k <= q/2
CenteredReduction := true;      // if true, use s1=0 reduction when p does not divide k

// -----------------------------------------------------------------------------
// Small utilities
// -----------------------------------------------------------------------------

function PairIndex(i, j, q)
    // i,j are 1-based element indices. The returned pair index is also 1-based.
    return (i - 1)*q + j;
end function;

function Popcount(n)
    // Number of 1-bits of a nonnegative integer n.
    c := 0;
    while n ne 0 do
        if BitwiseAnd(n, 1) ne 0 then
            c +:= 1;
        end if;
        n := ShiftRight(n, 1);
    end while;
    return c;
end function;

function FirstSetBit(n)
    // Return the 0-based position of the first 1-bit of n.
    // Assumes n > 0.
    pos := 0;
    while BitwiseAnd(n, 1) eq 0 do
        n := ShiftRight(n, 1);
        pos +:= 1;
    end while;
    return pos;
end function;

function FiniteFieldOrders(max_q, min_q, min_char)
    // Return triples <q,p,m> with q=p^m in [min_q,max_q] and p >= min_char.
    out := [];
    for p in [min_char..max_q] do
        if IsPrime(p) then
            q := p;
            m := 1;
            while q le max_q do
                if q ge min_q then
                    Append(~out, <q, p, m>);
                end if;
                m +:= 1;
                q *:= p;
            end while;
        end if;
    end for;
    Sort(~out);
    return out;
end function;

function BuildFieldData(q, p, m)
    /*
        Build GF(q), an ordered element list, and lookup tables.

        Element indexing convention:

        m = 1:
            index i represents F!(i-1), so bit position i-1 is the element i-1.

        m = 2:
            index a + p*b + 1 represents a + b*alpha, where alpha=F.1.
            This makes additive shifts easy: add delta=(da,db) by rotating rows.
    */

    F := GF(q);

    if m eq 1 then
        Elems := [ F!a : a in [0..p-1] ];
    else
        alpha := F.1;
        Elems := [ F!a + (F!b)*alpha : b in [0..p-1], a in [0..p-1] ];
    end if;

    IndexMap := AssociativeArray();
    for i in [1..q] do
        IndexMap[Elems[i]] := i;
    end for;

    AddIdx := [ [ IndexMap[Elems[i] + Elems[j]] : j in [1..q] ] : i in [1..q] ];
    SqIdx  := [ IndexMap[Elems[i]^2] : i in [1..q] ];
    CuIdx  := [ IndexMap[Elems[i]^3] : i in [1..q] ];

    // PairShift[xi][idx] sends (s1,s2) to (s1+x, s2+x^2).
    PairShift := [* *];
    for xi in [1..q] do
        x2i := SqIdx[xi];
        table := [ 0 : idx in [1..q*q] ];
        for s1i in [1..q] do
            ns1i := AddIdx[s1i][xi];
            for s2i in [1..q] do
                idx := PairIndex(s1i, s2i, q);
                ns2i := AddIdx[s2i][x2i];
                table[idx] := PairIndex(ns1i, ns2i, q);
            end for;
        end for;
        Append(~PairShift, table);
    end for;

    return F, Elems, AddIdx, SqIdx, CuIdx, PairShift;
end function;

function ShiftCubicMask(mask, delta_idx, p, m, q, full_mask, row_mask)
    /*
        Return the mask obtained from mask by adding delta to each cubic sum.

        The input delta is given by its 1-based element index delta_idx.
        Bit position c represents the field element with index c+1.

        For prime fields this is one cyclic rotation of q=p bits.
        For degree-2 fields, the additive group is GF(p)^2, so the mask is
        viewed as p rows of p bits and shifted by (da,db).
    */

    if mask eq 0 or delta_idx eq 1 or mask eq full_mask then
        return mask;
    end if;

    if m eq 1 then
        d := delta_idx - 1;
        return BitwiseAnd(BitwiseOr(ShiftLeft(mask, d), ShiftRight(mask, p - d)), full_mask);
    end if;

    // m = 2 case: index a + p*b + 1 represents a + b*alpha.
    delta := delta_idx - 1;
    da := delta mod p;
    db := delta div p;

    out := 0;
    for b in [0..p-1] do
        row := BitwiseAnd(ShiftRight(mask, p*b), row_mask);
        if row ne 0 then
            if da ne 0 then
                row := BitwiseAnd(BitwiseOr(ShiftLeft(row, da), ShiftRight(row, p - da)), row_mask);
            end if;
            out := BitwiseOr(out, ShiftLeft(row, p*((b + db) mod p)));
        end if;
    end for;

    return out;
end function;

function FirstMissingFull(DP, k, q, full_mask)
    /*
        Full check: all s1,s2,s3.
        Return:
            ok, first_s1_idx, first_s2_idx, first_s3_idx, total_missing
        The indices are 1-based element indices. They are meaningful only if ok=false.
    */

    total_missing := 0;
    have_first := false;
    fs1 := 1;
    fs2 := 1;
    fs3 := 1;

    layer := DP[k + 1];
    for idx in [1..q*q] do
        miss := BitwiseAnd(BitwiseNot(layer[idx]), full_mask);
        if miss ne 0 then
            if not have_first then
                pos := FirstSetBit(miss);
                fs1 := ((idx - 1) div q) + 1;
                fs2 := ((idx - 1) mod q) + 1;
                fs3 := pos + 1;
                have_first := true;
            end if;
            total_missing +:= Popcount(miss);
        end if;
    end for;

    return not have_first, fs1, fs2, fs3, total_missing;
end function;

function FirstMissingCentered(DP, k, q, full_mask)
    /*
        Centered check for p not dividing k.
        Only pairs (s1,s2)=(0,u) are tested.

        Return:
            ok, first_s1_idx, first_s2_idx, first_s3_idx,
            central_missing, implied_total_missing
    */

    central_missing := 0;
    have_first := false;
    fs1 := 1;
    fs2 := 1;
    fs3 := 1;

    layer := DP[k + 1];
    for uidx in [1..q] do
        idx := PairIndex(1, uidx, q);     // s1 = 0 has index 1
        miss := BitwiseAnd(BitwiseNot(layer[idx]), full_mask);
        if miss ne 0 then
            if not have_first then
                pos := FirstSetBit(miss);
                fs1 := 1;
                fs2 := uidx;
                fs3 := pos + 1;
                have_first := true;
            end if;
            central_missing +:= Popcount(miss);
        end if;
    end for;

    return not have_first, fs1, fs2, fs3, central_missing, central_missing*q;
end function;

function RunDP(q, p, m, max_k)
    /*
        Build the DP table.

        DP[t+1][PairIndex(s1,s2,q)] is a q-bit integer mask.
        The bit for s3 is on iff a t-subset with moments (s1,s2,s3) exists.
    */

    F, Elems, AddIdx, SqIdx, CuIdx, PairShift := BuildFieldData(q, p, m);

    full_mask := 2^q - 1;
    row_mask := 2^p - 1;

    DP := [ [ 0 : idx in [1..q*q] ] : t in [0..max_k] ];
    DP[1][PairIndex(1, 1, q)] := 1;       // empty subset has (0,0,0)

    TableFull := [ false : t in [0..max_k] ];

    q2 := q*q;

    for seen in [1..q] do
        xi := seen;
        x3i := CuIdx[xi];
        upper := Minimum(seen, max_k);
        shift_pair := PairShift[xi];

        // Descend in t so that xi is used at most once.
        for t := upper to 1 by -1 do
            if not TableFull[t + 1] then
                for idx in [1..q2] do
                    mask := DP[t][idx];
                    if mask ne 0 then
                        j := shift_pair[idx];
                        if DP[t + 1][j] ne full_mask then
                            shifted := ShiftCubicMask(mask, x3i, p, m, q, full_mask, row_mask);
                            DP[t + 1][j] := BitwiseOr(DP[t + 1][j], shifted);
                        end if;
                    end if;
                end for;

                // If all fibers are already full, future updates are unnecessary.
                is_full := true;
                for idx in [1..q2] do
                    if DP[t + 1][idx] ne full_mask then
                        is_full := false;
                        break;
                    end if;
                end for;
                if is_full then
                    TableFull[t + 1] := true;
                end if;
            end if;
        end for;
    end for;

    return DP, Elems, full_mask;
end function;

function CheckField(q, p, m, min_k, use_symmetry, centered_reduction)
    max_k := q - min_k;
    if use_symmetry then
        max_k := q div 2;
    end if;

    t0 := Cputime();
    DP, Elems, full_mask := RunDP(q, p, m, max_k);
    secs := Cputime(t0);

    field_ok := true;
    bad_lines := [* *];

    for k in [min_k..max_k] do
        can_center := centered_reduction and ((k mod p) ne 0);

        if can_center then
            ok, s1i, s2i, s3i, central_missing, implied_total_missing :=
                FirstMissingCentered(DP, k, q, full_mask);
            if not ok then
                field_ok := false;
                Append(~bad_lines,
                    <k, true, Elems[s1i], Elems[s2i], Elems[s3i],
                     central_missing, implied_total_missing>);
            end if;
        else
            ok, s1i, s2i, s3i, total_missing := FirstMissingFull(DP, k, q, full_mask);
            if not ok then
                field_ok := false;
                Append(~bad_lines,
                    <k, false, Elems[s1i], Elems[s2i], Elems[s3i],
                     total_missing, total_missing>);
            end if;
        end if;
    end for;

    return field_ok, bad_lines, secs, max_k;
end function;

// -----------------------------------------------------------------------------
// Main program
// -----------------------------------------------------------------------------

fields := FiniteFieldOrders(MaxQ, MinQ, Maximum(MinChar, 5));

printf "Fields to check: ";
for i in [1..#fields] do
    triple := fields[i];
    printf "%o=%o^%o", triple[1], triple[2], triple[3];
    if i lt #fields then
        printf ", ";
    end if;
end for;
printf "\n\n";

all_ok := true;
total_secs := 0.0;

for triple in fields do
    q := triple[1];
    p := triple[2];
    m := triple[3];

    ok, bad_lines, secs, max_k := CheckField(q, p, m, MinK, UseSymmetry, CenteredReduction);
    total_secs +:= secs;

    if m eq 1 then
        printf "q=%o, GF(%o), k=%o..%o, DP time=%o seconds\n", q, q, MinK, max_k, secs;
    else
        printf "q=%o, GF(%o^%o), k=%o..%o, DP time=%o seconds\n", q, p, m, MinK, max_k, secs;
    end if;

    if ok then
        printf "  OK\n\n";
    else
        all_ok := false;
        printf "  FAIL\n";
        for line in bad_lines do
            k := line[1];
            centered := line[2];
            b1 := line[3];
            b2 := line[4];
            b3 := line[5];
            count1 := line[6];
            count2 := line[7];

            if centered then
                printf "    k=%o: centered s1=0 check; first missing centered triple=(%o,%o,%o); central missing=%o; implied total missing=%o\n",
                    k, b1, b2, b3, count1, count2;
            else
                printf "    k=%o: full check; first missing triple=(%o,%o,%o); total missing=%o\n",
                    k, b1, b2, b3, count1;
            end if;
        end for;
        printf "\n";
    end if;
end for;

if all_ok then
    printf "OVERALL RESULT: all checked cases passed. Total DP time=%o seconds\n", total_secs;
else
    printf "OVERALL RESULT: at least one checked case failed. Total DP time=%o seconds\n", total_secs;
end if;
