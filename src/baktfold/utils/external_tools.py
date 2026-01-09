"""
Originally taken from Michael Hall's tbpore https://github.com/mbhall88/tbpore/blob/main/tbpore/external_tools.py

Also used by a variety of other tools (Dnaapler, Plassembler, Pharokka)

"""

import hashlib
import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import click
from loguru import logger
 

class ExternalTool:
    """
    Class for running external tools.

    Args:
      tool (str): The path to the tool to run.
      input (str): The input file.
      output (str): The output file.
      params (str): The parameters to pass to the tool.
      logdir (Path): The directory to store log files.

    Attributes:
      command (List[str]): The command to run.
      out_log (str): The path to the stdout log file.
      err_log (str): The path to the stderr log file.

    Examples:
      >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
      >>> tool.command
      ["tool", "params", "output", "input"]
      >>> tool.out_log
      "logdir/tool_1234567890abcdef1234567890abcdef.out"
      >>> tool.err_log
      "logdir/tool_1234567890abcdef1234567890abcdef.err"
    """
    def __init__(self, tool: str, input: str, output: str, params: str, logdir: Path):
        """
        Initializes an ExternalTool object.

        Args:
          tool (str): The path to the tool to run.
          input (str): The input file.
          output (str): The output file.
          params (str): The parameters to pass to the tool.
          logdir (Path): The directory to store log files.

        Attributes:
          command (List[str]): The command to run.
          out_log (str): The path to the stdout log file.
          err_log (str): The path to the stderr log file.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool.command
          ["tool", "params", "output", "input"]
          >>> tool.out_log
          "logdir/tool_1234567890abcdef1234567890abcdef.out"
          >>> tool.err_log
          "logdir/tool_1234567890abcdef1234567890abcdef.err"
        """
        logdir = Path(logdir)   
        self.command: List[str] = self._build_command(tool, input, output, params)
        Path(logdir).mkdir(parents=True, exist_ok=True)
        command_hash = hashlib.sha256(self.command_as_str.encode("utf-8")).hexdigest()
        tool_name = Path(tool).name
        logfile_prefix: Path = logdir / f"{tool_name}_{command_hash}"
        self.out_log = f"{logfile_prefix}.out"
        self.err_log = f"{logfile_prefix}.err"

    @property
    def command_as_str(self) -> str:
        """
        Returns the command as a string.

        Returns:
          str: The command as a string.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool.command_as_str
          "tool params output input"
        """
        return shlex.join(self.command)

    @staticmethod
    def _build_command(tool: str, input: str, output: str, params: str) -> List[str]:
        """
        Builds the command to run.

        Args:
          tool (str): The path to the tool to run.
          input (str): The input file.
          output (str): The output file.
          params (str): The parameters to pass to the tool.

        Returns:
          List[str]: The command to run.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool._build_command("tool", "input", "output", "params")
          ["tool", "params", "output", "input"]
        """
        # note: shlex.join does not allow us to shlex.split() later
        # this is explicitly a " ".join()
        command = " ".join([tool, params, output, input])
        escaped_command = shlex.split(command)
        return escaped_command

    def run(self) -> None:
        """
        Runs the tool.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool.run()
        """
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            print(f"Command line: {self.command_as_str}", file=stderr_fh)
            logger.info(f"Started running {self.command_as_str} ...")
            self._run_core(self.command, stdout_fh=stdout_fh, stderr_fh=stderr_fh)
            logger.info(f"Done running {self.command_as_str}")

    """
    stream to terminal (aria2c) so the user knows how long it is taking
    """

    def run_stream(self) -> None:
        """
        Runs the tool and streams the output to the terminal.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool.run_stream()
        """
        with open(self.out_log, "w") as stdout_fh, open(self.err_log, "w") as stderr_fh:
            print(f"Command line: {self.command_as_str}", file=stderr_fh)
            logger.info(f"Started running {self.command_as_str} ...")

            process = subprocess.Popen(
                self.command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                bufsize=1,
                universal_newlines=True,
            )

            for line in process.stdout:
                print(line, end="")         # Live output to terminal
                stdout_fh.write(line)       # Also write to stdout log

            process.stdout.close()
            return_code = process.wait()

            logger.info(f"Done running {self.command_as_str}")

            if return_code != 0:
                raise subprocess.CalledProcessError(return_code, self.command)


    @staticmethod
    def _run_core(command: List[str], stdout_fh, stderr_fh) -> None:
        """
        Runs the tool.

        Args:
          command (List[str]): The command to run.
          stdout_fh: The file handle to write stdout to.
          stderr_fh: The file handle to write stderr to.

        Examples:
          >>> tool = ExternalTool("tool", "input", "output", "params", "logdir")
          >>> tool._run_core(["tool", "params", "output", "input"], stdout_fh, stderr_fh)
        """
        subprocess.check_call(command, stdout=stdout_fh, stderr=stderr_fh)

    @staticmethod
    def run_tools(
        tools_to_run: Tuple["ExternalTool", ...], ctx: Optional[click.Context] = None
    ) -> None:
        """
        Runs a list of tools.

        Args:
          tools_to_run (Tuple[ExternalTool]): The list of tools to run.
          ctx (Optional[click.Context]): The click context.

        Examples:
          >>> tool1 = ExternalTool("tool1", "input1", "output1", "params1", "logdir")
          >>> tool2 = ExternalTool("tool2", "input2", "output2", "params2", "logdir")
          >>> ExternalTool.run_tools((tool1, tool2))
          >>> ExternalTool.run_tools((tool1, tool2), ctx)
        """
        for tool in tools_to_run:
            try:
                tool.run()
            except subprocess.CalledProcessError as error:
                logger.error(
                    f"Error calling {tool.command_as_str} (return code {error.returncode})"
                )
                logger.error(f"Please check stdout log file: {tool.out_log}")
                logger.error(f"Please check stderr log file: {tool.err_log}")
                logger.error("Temporary files are preserved for debugging")
                logger.error("Exiting...")

                if ctx:
                    ctx.exit(1)
                else:
                    sys.exit(1)

    """
    Only one toolf
    """

    @staticmethod
    def run_tool(tool: "ExternalTool", ctx: Optional[click.Context] = None) -> None:
        """
        Runs the given external tool.

        Args:
          tool (ExternalTool): The external tool to run.
          ctx (Optional[click.Context]): The click context to use. Defaults to None.

        Returns:
          None.

        Raises:
          subprocess.CalledProcessError: If there is an error calling the external tool.

        Examples:
          >>> tool = ExternalTool()
          >>> ExternalTool.run_tool(tool)
          None
        """
        try:
            tool.run()
        except subprocess.CalledProcessError as error:
            logger.error(
                f"Error calling {tool.command_as_str} (return code {error.returncode})"
            )
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error("Temporary files are preserved for debugging")
            logger.error("Exiting...")

            if ctx:
                ctx.exit(1)
            else:
                sys.exit(1)


    """
    Only download - so can print the aria2c output to screen
    """

    @staticmethod
    def run_download(tool: "ExternalTool", ctx: Optional[click.Context] = None) -> None:
        """
        Runs the given external tool and prints the aria2c output to the screen.

        Args:
          tool (ExternalTool): The external tool to run.
          ctx (Optional[click.Context]): The click context to use. Defaults to None.

        Returns:
          None.

        Raises:
          subprocess.CalledProcessError: If there is an error calling the external tool.

        Examples:
          >>> tool = ExternalTool()
          >>> ExternalTool.run_download(tool)
          None
        """
        try:
            tool.run_stream()
        except subprocess.CalledProcessError as error:
            logger.error(
                f"Error calling {tool.command_as_str} (return code {error.returncode})"
            )
            logger.error(f"Please check stdout log file: {tool.out_log}")
            logger.error(f"Please check stderr log file: {tool.err_log}")
            logger.error("Temporary files are preserved for debugging")
            logger.error("Exiting...")

            if ctx:
                ctx.exit(1)
            else:
                sys.exit(1)
