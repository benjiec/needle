LLM, please following these rules

1. When answering a question, don't assume you need to make changes
immediately. Propose solution, get confirmation, then make the change.

2. Don't fix the a bug just to get rid of the error, but remember the purpose
of the changes made that may have caused the bug and preserve that purpose.

3. Prefer concise answers to simple questions, rather than expanding on an
answer with possibilities and providing lots of text to read.

4. Please add new unit tests whenever you are adding new features. New tests
should follow Python unittest module convention (e.g. subclass from
unittest.TestCase). IMPORTANT - place all tests for a .py file into the same
test .py file, i.e. don't put different behavior tests for the same .py file
into multiple test files. Run unit tests using `scripts/run-unit-tests`.

5. Don't leave comments that only pertain to an action you took, e.g. to undo
something. Comments should describe the current state of the code or
configuration.

6. Follow DRY model in code development. Abstract common code into functions
that can be re-called as much as possible.

7. DO NOT seed the random number generator in code. Put random seeding in unit
tests instead to get consistent behavior.

9. When removing old code we are not using, which you should almost always do,
please remove all the code including the function signature/definitino, and
don't just leave an empty "pass" or "return None" as a placeholder. And remove
tests for logic we don't use anymore. And, DO NOT leave a comment saying you
removed something.

10. Don't worry about backwards compatibility. Propose the most
concise/clean/DRY solution and then propose possible backwards compatibility
approaches. But don't just create code and leave code around for backward
compatibility reasons.

11. Do not use try/except unless catching a specific, documented exception
type. Never use bare except or catch broad types like Exception/BaseException.
For example, please DO NOT use "except:" or "except Exception:".

12. Do not use single line functions that just returns outcome of another
function with the exact same arguments. Just go and change all the use of
former to the latter.

13. Do not leave empty functions around if they are not used anymore. E.g. "def
hello(): pass" or "def hello(): return". Just remove code we don't use.

14. Do not use this style "getattr(self, <field name>, None)" - just set the
field in the class and use self.<field name>.
