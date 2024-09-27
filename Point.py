
class point:

    def __init__(self, C_Si, C_Mo, C_N, C_O2) -> None:
        self.C_Si = C_Si
        self.C_Mo = C_Mo
        self.C_N = C_N
        self.C_O2 = C_O2

    @staticmethod
    def loss(p1, p2):
        err = abs(p1.C_Si - p2.C_Si) + abs(p1.C_Mo - p2.C_Mo) + abs(p1.C_N - p2.C_N) + abs(p1.C_O2 - p2.C_O2)
        return err

    def __repr__(self) -> str:
        return f"{self.C_Si, self.C_Mo, self.C_N, self.C_O2}, "

    def __str__(self) -> str:
        return f"Si: {self.C_Si}\nMo: {self.C_Mo}\nN: {self.C_N}\nO: {self.C_O2}"
