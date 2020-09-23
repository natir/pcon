import pcon

def test_count():
    a = pcon.Counter(5)
    
    a.count_fasta("../data/test.fasta", 10)

    assert a.get(108) == 18

    a.serialize("python_counter.pcon", 0)

    b = pcon.Counter.deserialize("python_counter.pcon")

    assert a.get(108) == b.get(108)

    b.inc(108)

    assert b.get(108) == 19

def test_csv():
    a = pcon.Counter(5)

    pcon.dump.csv(a, 0, "python_counter.csv")

def test_solid():
    a = pcon.Counter(5)

    pcon.dump.solid(a, 0, "python_counter.solid")

def test_spectrum():
    a = pcon.Counter(5)

    pcon.dump.spectrum(a, "python_counter.spectrum.csv")
