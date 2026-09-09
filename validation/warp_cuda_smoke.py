import warp as wp


@wp.kernel
def add_one(x: wp.array(dtype=float)):
    i = wp.tid()
    x[i] += 1.0


def main() -> None:
    wp.init()

    x = wp.array([1.0, 2.0, 3.0], dtype=float, device="cuda:0")
    wp.launch(add_one, dim=3, inputs=[x], device="cuda:0")
    wp.synchronize()

    result = x.numpy()
    print(result)

    if result.tolist() != [2.0, 3.0, 4.0]:
        raise RuntimeError(f"unexpected result: {result}")


if __name__ == "__main__":
    main()
