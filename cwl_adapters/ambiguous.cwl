cwlVersion: v1.2

class: CommandLineTool

requirements:
  InlineJavascriptRequirement: {}

baseCommand: "true"

# This CommandLineTool is solely to demonstrate the --inference_use_naming_conventions feature.
# Each output simply returns the corresponding input.
# Since both outputs are of type string, inference (by default) cannot
# distinguish between the two outputs, and will always choose uniquename2
# since the outputs are considered ordered and matched from botton to top.

inputs:
  uniquename1:
    type: string
    format: "someformat"

  uniquename2:
    type: string
    format: "someformat"

outputs:
  uniquename1:
    type: string
    # outputs of type: string cannot have formats
    outputBinding:
      outputEval: $(inputs.uniquename1)

  uniquename2:
    type: string
    # outputs of type: string cannot have formats
    outputBinding:
      outputEval: $(inputs.uniquename2)