# MSI
## Current Project Structure
``` Bash
MSI
├── cti_core
│   ├── cti_combine.py
│   ├── cti_processor.py
│   ├── soln2cti.py
│   └── testind
├── data
├── instruments
│   └── shock_tube.py
├── README.md
├── simulations
│   ├── free_flame.py
│   └── simulations.py
├── tests
│   
└── utilities
    └── cti_remover.py
```

### Usage:

Example: CTI pruning in a python prompt-
``` Python
Python 3.6.5 |Anaconda, Inc.| (default, Apr 29 2018, 16:14:56) 
[GCC 7.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import MSI.cti_core.cti_processor as pr
>>> test_processor = pr.Processor('./MSI/data/test_data/FFCM1.cti')
>>> test_processor.prune('./MSI/cti_core/testind')
Error on index 5say
 Skipping index
remove index 1, reaction H + O2 <=> O + OH
remove index 2, reaction H2 + O <=> H + OH
remove index 3, reaction H2 + O <=> H + OH
>>> test_processor.prune(3)
remove index 3, reaction H2 + M <=> 2 H + M
>>> test_processor.prune([1,4,77])
remove index 1, reaction H + O2 <=> O + OH
remove index 4, reaction H2 + HE <=> 2 H + HE
remove index 77, reaction CH2(S) + H2O2 <=> CH3O + OH
>>> test_processor.write_soln_to_file()
'./MSI/data/test_data/FFCM1_processed.cti'
>>> exit()
```

Example: active parameter setting, editing, and getting in a python prompt-
``` Python
Python 3.6.5 |Anaconda, Inc.| (default, Apr 29 2018, 16:14:56) 
[GCC 7.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from MSI.cti_core import cti_processor
>>> test = cti_processor.Processor('MSI/data/test_data/FFCM1.cti')
>>> test.set_default_parameters()
>>> test.get_active_parameter(1)
<MSI.cti_core.cti_processor.active_parameter object at 0x7fc90cf21518>
>>> test.get_active_parameter(1,1)
Reaction H + O2 <=> O + OH:
Type: ElementaryReaction
dels: [0.0, 0.0, 0.0]
h_dels: []
l_dels: []
 rate_list: []
<MSI.cti_core.cti_processor.active_parameter object at 0x7fc90cf21518>
>>> test.add_active_parameter(1,'ElementaryReaction',dels=[.6,.6,.6])
True
>>> test.get_active_parameter(1,1)
Reaction H + O2 <=> O + OH:
Type: ElementaryReaction
dels: [0.6, 0.6, 0.6]
h_dels: []
l_dels: []
 rate_list: []
<MSI.cti_core.cti_processor.active_parameter object at 0x7fc8ff7cb978>
>>> exit()
```
Example: active parameter setting with processor initialization in a python prompt
``` Python
Python 3.6.5 |Anaconda, Inc.| (default, Apr 29 2018, 16:14:56) 
[GCC 7.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from MSI.cti_core import cti_processor
>>> test = cti_processor.Processor('MSI/data/test_data/FFCM1.cti')
>>> test.get_active_parameter(1)
Error: active parameter dictionary empty. Please initialize or add a parameter to the dictionary
-1
>>> test = cti_processor.Processor('MSI/data/test_data/FFCM1.cti',1)
>>> test.get_active_parameter(1)
<MSI.cti_core.cti_processor.active_parameter object at 0x7f04fe6cbc18>
>>> test.get_active_parameter(1,1)
Reaction H + O2 <=> O + OH:
Type: ElementaryReaction
dels: [0.0, 0.0, 0.0]
h_dels: []
l_dels: []
 rate_list: []
<MSI.cti_core.cti_processor.active_parameter object at 0x7f04fe6cbc18>
>>> exit()
```

Example: writing active parameters to file in a python prompt:
``` Python
Python 3.6.5 |Anaconda, Inc.| (default, Apr 29 2018, 16:14:56) 
[GCC 7.2.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from MSI.cti_core import cti_processor
>>> test = cti_processor.Processor('MSI/data/test_data/FFCM1.cti',1)
>>> test.add_active_parameter(1,'ElementaryReaction',dels=[.6,.6,.6])
True
>>> test.write_parameters()
'MSI/data/test_data/FFCM1_processed.param'
>>> exit()
```
Data is formatted like this:
```
Type: ThreeBodyReaction
dels: [0.0, 0.0, 0.0]
Reaction 14: '2 H2O <=> H + H2O + OH'
Type: ElementaryReaction
dels: [0.0, 0.0, 0.0]
Reaction 15: 'H + O2 (+M) <=> HO2 (+M)'
Type: FalloffReaction
h_dels: [0.0, 0.0, 0.0]
l_dels: [0.0, 0.0, 0.0]
```
