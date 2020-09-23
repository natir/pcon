import pcon

def test_solid():
    a = pcon.Counter(5)
    a.deserialize("python_counter.pcon")

    b = pcon.Solid.from_counter(a, 20)

    assert b.get(108) == False

    b.serialize("python_solid.pcon")

    c = pcon.Solid.deserialize("python_solid.pcon")

    assert c.get(108) == False

    c.set(108, True)

    assert c.get(108) == True
