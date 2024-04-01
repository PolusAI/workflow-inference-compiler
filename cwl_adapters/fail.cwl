cwlVersion: v1.0

class: CommandLineTool

label: Randomly fail with a given probability. Useful for testing low probability failures when scattering.

doc: |
  Randomly fail with a given probability. Useful for testing low probability failures when scattering.

# NOTE: no explicit baseCommand; return true or false in javascript.
# See https://www.commonwl.org/user_guide/topics/expressions.html
arguments:
- valueFrom: |
    ${
      if (Math.random() < inputs.with_probability) {
        return "false"; // NOTE: These must be strings!
      } else {
        return "true"; // NOTE: These must be strings!
      }
    }
# I do not recommend trying to do this in bash, because between bash itself and
# CWL, the string escaping is nearly impossible to get right.

requirements:
  InlineJavascriptRequirement: {}

inputs:
  with_probability:
    type: float
    default: 0.01  # i.e. 1%

outputs:
  failout:
    type: string
    outputBinding:
      outputEval: $("not fail")
