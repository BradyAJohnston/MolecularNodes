# Docstring Checklist

Run through this checklist for every public symbol in the package.

## Functions and methods

- [ ] One-line summary in imperative mood ("Return the...", not
      "Returns the...")
- [ ] Extended summary if behavior is non-obvious
- [ ] Every parameter documented with description
- [ ] Type annotations present in the signature (not just docstring)
- [ ] `Returns` section with type and description
- [ ] `Raises` section for each exception that can be raised
- [ ] At least one `Examples` entry — prefer Quarto `{python}` cells
      with prose; `>>>` prompts are acceptable for simple cases
- [ ] `%seealso` directive for closely related functions
- [ ] No private parameters documented (leading `_`)
- [ ] Default values mentioned in parameter descriptions

## Classes

- [ ] Class-level docstring (not on `__init__`)
- [ ] One-line summary describing what the class represents
- [ ] Constructor parameters documented under `Parameters`
- [ ] Key public methods mentioned in extended summary or `Notes`
- [ ] `Examples` showing basic instantiation and usage

## Properties

- [ ] One-line summary describing what the property returns
- [ ] Type annotation on the `@property` return
- [ ] `Returns` section if the value is complex

## Module-level constants

- [ ] Docstring or inline comment explaining the constant's purpose
- [ ] Type annotation

## Common issues to flag

- Docstring says "Returns" but function returns `None`
- Parameter name in docstring doesn't match signature
- Missing blank line between summary and extended summary
- Inconsistent style (NumPy in some files, Google in others)
- Examples that import from internal modules (`_private`)
- Stale parameter documentation after a signature change
