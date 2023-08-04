import logging
import re

from src.covsonar.logging import LoggingConfigurator


def test_logging(logger, caplog, capsys):
    logger_configurator = LoggingConfigurator()
    logger_configurator.set_debug_mode(False)

    # Test non-debug mode
    with caplog.at_level(logging.INFO, logger=logger.name):
        # Debug message (ignored in non-debug mode)
        logger.debug("This is a debug message (non-debug mode).")
        assert len(caplog.records) == 0

        # Info message
        logger.info("This is an info message (non-debug mode).")
        assert "This is an info message (non-debug mode)." == caplog.records[-1].message
        captured = capsys.readouterr()
        assert "This is an info message (non-debug mode)." in captured.out

        # Error message
        logger.error("This is an error message (non-debug mode).")
        assert (
            "This is an error message (non-debug mode)." == caplog.records[-1].message
        )
        captured = capsys.readouterr()
        assert "ERROR: This is an error message (non-debug mode)." in captured.err

        # Critical message
        logger.critical("This is a critical message (non-debug mode).")
        assert (
            "This is a critical message (non-debug mode)." == caplog.records[-1].message
        )
        captured = capsys.readouterr()
        assert "CRITICAL: This is a critical message (non-debug mode)." in captured.err


def test_logging_debug(logger, caplog, capsys):
    logger_configurator = LoggingConfigurator()
    logger_configurator.set_debug_mode(True)

    # Define a common pattern for timestamp, logger name, file name, and line number
    common_pattern = (
        r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2},\d{3} - "
        + rf"{logger.name} - [\w_]+.py:\d+ - "
    )

    with caplog.at_level(logging.DEBUG, logger=logger.name):
        # Debug message
        logger.debug("This is a debug message (debug mode).")
        pattern = re.compile(
            common_pattern + r"DEBUG: This is a debug message \(debug mode\)\."
        )
        assert pattern.search(capsys.readouterr().err)

        # Info message
        logger.info("This is an info message (debug mode).")
        pattern = re.compile(
            common_pattern + r"INFO: This is an info message \(debug mode\)\."
        )
        assert pattern.search(capsys.readouterr().out)

        # Error message
        logger.error("This is an error message (debug mode).")
        pattern = re.compile(
            common_pattern + r"ERROR: This is an error message \(debug mode\)\."
        )
        assert pattern.search(capsys.readouterr().err)

        # Critical message
        logger.critical("This is a critical message (debug mode).")
        pattern = re.compile(
            common_pattern + r"CRITICAL: This is a critical message \(debug mode\)\."
        )
        assert pattern.search(capsys.readouterr().err)
