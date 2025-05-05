def pprint_format_as_set(
    values: list,
    keys: list[str] = None
) -> None:
    """
    Pretty prints a set of values and keys in mathematical set format.

    Args:
        values (list):
            Values to print.
        keys (list[str]):
            Keys corresponding with values to print. Defaults to None. 
    """

    if values:
        output = f'{{{tuple(keys)}:'.replace("'", "") if keys else '{'
        for value in values:
            output += f'\n {tuple(value)}'
        output += '\n}'
    else:
        output = '{}'
    print(output)
