from Utilities.ElectricalObjects import NodeFactory, LineFactory, Grid


def main():
    # Creazione dei nodi
    slack = NodeFactory.create_slack_bus(voltage=1.05, name="Slack")
    gen1 = NodeFactory.create_generator_bus(
        voltage=1.05, p_gen=0.5, q_min=-0.5, q_max=1.0, name="Gen1"
    )
    gen2 = NodeFactory.create_generator_bus(
        voltage=1.07, p_gen=0.6, q_min=-0.5, q_max=1.5, name="Gen2"
    )
    load1 = NodeFactory.create_load_bus(p_load=0.7, q_load=0.7, name="Load1")
    load2 = NodeFactory.create_load_bus(p_load=0.7, q_load=0.7, name="Load2")
    load3 = NodeFactory.create_load_bus(p_load=0.7, q_load=0.7, name="Load3")

    # Creazione delle linee
    line1 = LineFactory.create_transmission_line(
        slack, gen1, r=0.1, x=0.2, b=0.04, name="Line1"
    )
    line2 = LineFactory.create_transmission_line(
        slack, load1, r=0.05, x=0.2, b=0.04, name="Line2"
    )
    line3 = LineFactory.create_transmission_line(
        slack, load2, r=0.08, x=0.3, b=0.06, name="Line3"
    )
    line4 = LineFactory.create_transmission_line(
        gen1, gen2, r=0.05, x=0.25, b=0.06, name="Line4"
    )
    line5 = LineFactory.create_transmission_line(
        gen1, load1, r=0.05, x=0.1, b=0.02, name="Line5"
    )
    line6 = LineFactory.create_transmission_line(
        gen1, load2, r=0.1, x=0.3, b=0.04, name="Line6"
    )
    line7 = LineFactory.create_transmission_line(
        gen1, load3, r=0.07, x=0.2, b=0.05, name="Line7"
    )
    line8 = LineFactory.create_transmission_line(
        gen2, load2, r=0.12, x=0.26, b=0.05, name="Line8"
    )
    line9 = LineFactory.create_transmission_line(
        gen2, load3, r=0.02, x=0.1, b=0.02, name="Line9"
    )
    line10 = LineFactory.create_transmission_line(
        load1, load2, r=0.2, x=0.4, b=0.08, name="Line10"
    )
    line11 = LineFactory.create_transmission_line(
        load2, load3, r=0.1, x=0.3, b=0.06, name="Line11"
    )

    # Creazione della rete
    nodes = [slack, gen1, gen2, load1, load2, load3]
    lines = [
        line1,
        line2,
        line3,
        line4,
        line5,
        line6,
        line7,
        line8,
        line9,
        line10,
        line11,
    ]

    grid = Grid(nodes=nodes, lines=lines)
    grid.loadflow()


if __name__ == "__main__":
    main()
