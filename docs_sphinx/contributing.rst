Contributing to ClearWater Modules
==================================

We welcome contributions to the ClearWater Modules project! This guide will help you get started with contributing to the codebase.

Getting Started
---------------

Development Setup
~~~~~~~~~~~~~~~~~

1. Fork the repository on GitHub
2. Clone your fork locally:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/ClearWater-modules.git
      cd ClearWater-modules

3. Create a development environment:

   .. code-block:: bash

      conda env create -f environment.yml
      conda activate clearwater-modules

4. Install in development mode:

   .. code-block:: bash

      conda develop .

5. Create a feature branch:

   .. code-block:: bash

      git checkout -b feature-name

Code Standards
--------------

Python Style
~~~~~~~~~~~~

We follow PEP 8 style guidelines with the following specifications:

- Line length: 88 characters (Black formatter default)
- Use type hints where appropriate
- Docstrings: NumPy style

Example:

.. code-block:: python

   def calculate_rate(
       concentration: float,
       temperature: float,
       rate_constant: float = 0.1
   ) -> float:
       """
       Calculate reaction rate using Arrhenius equation.
       
       Parameters
       ----------
       concentration : float
           Reactant concentration in mg/L
       temperature : float
           Water temperature in degrees Celsius
       rate_constant : float, optional
           Base rate constant at 20°C, by default 0.1
       
       Returns
       -------
       float
           Reaction rate in mg/L/day
       
       Notes
       -----
       Uses a temperature coefficient of 1.047
       """
       theta = 1.047
       rate = concentration * rate_constant * theta ** (temperature - 20)
       return rate

Code Organization
~~~~~~~~~~~~~~~~~

- Each module should have its own directory
- Separate concerns into different files:
  - ``constants.py``: Physical constants
  - ``state_variables.py``: Time-varying state
  - ``static_variables.py``: Parameters
  - ``processes.py``: Scientific calculations
  - ``model.py``: Main interface

Testing
-------

Writing Tests
~~~~~~~~~~~~~

All new features should include tests. We use pytest for testing:

.. code-block:: python

   # test_my_feature.py
   import pytest
   from clearwater_modules.my_module import my_function
   
   def test_my_function_basic():
       """Test basic functionality."""
       result = my_function(10, 20)
       assert result == 30
   
   def test_my_function_edge_case():
       """Test edge cases."""
       with pytest.raises(ValueError):
           my_function(-1, 20)
   
   @pytest.mark.parametrize("input,expected", [
       (0, 0),
       (1, 1),
       (10, 100),
   ])
   def test_my_function_parametrized(input, expected):
       """Test with multiple inputs."""
       result = my_function(input, input)
       assert result == expected

Running Tests
~~~~~~~~~~~~~

Run the test suite:

.. code-block:: bash

   # Run all tests
   pytest
   
   # Run specific test file
   pytest tests/test_my_feature.py
   
   # Run with coverage
   pytest --cov=clearwater_modules

Documentation
-------------

Docstrings
~~~~~~~~~~

All public functions and classes need docstrings:

.. code-block:: python

   class WaterQualityModel:
       """
       Base class for water quality models.
       
       This class provides the framework for implementing water quality
       simulation modules with consistent interfaces.
       
       Attributes
       ----------
       timestep : float
           Simulation timestep in seconds
       state : dict
           Current state variables
       
       Methods
       -------
       step(dt)
           Advance the model by one timestep
       get_state()
           Return current state dictionary
       """

Adding to Sphinx Docs
~~~~~~~~~~~~~~~~~~~~~

When adding new modules:

1. Create a new `.rst` file in `docs_sphinx/modules/`
2. Update the table of contents in `docs_sphinx/modules/index.rst`
3. Add API documentation in `docs_sphinx/api/`

Contribution Process
--------------------

Making Changes
~~~~~~~~~~~~~~

1. Make your changes in a feature branch
2. Add or update tests as needed
3. Update documentation
4. Run tests to ensure everything passes
5. Commit your changes with clear messages:

   .. code-block:: bash

      git add .
      git commit -m "Add feature: brief description
      
      Longer explanation of what changed and why"

Submitting Pull Requests
~~~~~~~~~~~~~~~~~~~~~~~~

1. Push your branch to your fork:

   .. code-block:: bash

      git push origin feature-name

2. Open a pull request on GitHub
3. Fill out the PR template with:
   
   - Description of changes
   - Related issues
   - Test results
   - Documentation updates

4. Address review comments
5. Once approved, we'll merge your contribution!

Code Review Process
~~~~~~~~~~~~~~~~~~~

All contributions go through code review. Reviewers will check for:

- Code quality and style
- Test coverage
- Documentation completeness
- Scientific accuracy
- Performance considerations

Types of Contributions
----------------------

Bug Reports
~~~~~~~~~~~

Report bugs at https://github.com/EcohydrologyTeam/ClearWater-modules/issues

Include:

- Your operating system
- Python version
- Detailed steps to reproduce
- Expected vs actual behavior
- Error messages and tracebacks

Feature Requests
~~~~~~~~~~~~~~~~

Submit feature requests as GitHub issues. Describe:

- The problem you're trying to solve
- How the feature would work
- Why it belongs in ClearWater

Documentation
~~~~~~~~~~~~~

Help improve documentation by:

- Fixing typos and clarifying text
- Adding examples
- Creating tutorials
- Improving API documentation

Scientific Contributions
~~~~~~~~~~~~~~~~~~~~~~~~

We especially welcome:

- New process algorithms
- Improved numerical methods
- Validation studies
- Case study examples

Development Guidelines
----------------------

Module Architecture
~~~~~~~~~~~~~~~~~~~

New modules should follow the established pattern:

.. code-block:: python

   # my_module/model.py
   from clearwater_modules.base import WaterQualityModule
   
   class MyModule(WaterQualityModule):
       def __init__(self):
           super().__init__()
           self.name = "MyModule"
       
       def initialize(self, state_vars, static_vars):
           """Initialize the module."""
           pass
       
       def step(self, timestep):
           """Advance one timestep."""
           pass

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~

- Use NumPy for array operations
- Avoid unnecessary loops
- Profile code for bottlenecks
- Consider memory usage for large simulations

Version Control
~~~~~~~~~~~~~~~

- Keep commits focused and atomic
- Write clear commit messages
- Reference issues in commits
- Keep feature branches up to date with main

Community
---------

Communication Channels
~~~~~~~~~~~~~~~~~~~~~~

- GitHub Issues: Bug reports and feature requests
- GitHub Discussions: General questions and discussions
- Email: clearwater@erdc.usace.army.mil

Code of Conduct
~~~~~~~~~~~~~~~

We are committed to providing a welcoming and inclusive environment. Please:

- Be respectful and considerate
- Welcome newcomers
- Focus on constructive criticism
- Respect differing viewpoints

Recognition
-----------

Contributors will be recognized in:

- The project README
- Release notes
- Scientific publications (for substantial contributions)

Thank you for contributing to ClearWater Modules!