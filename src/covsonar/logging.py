import configparser
import logging
import sys
from typing import Dict, Optional


class LoggingConfigurator:
    """A class to configure logging.

    INFO messages are directed to stdout, while DEBUG, WARNING, ERROR, and CRITICAL messages are directed to stderr.
    The format can be either simple or detailed based on the debug flag. Without debug, the level name is omitted for INFO messages.
    """

    class InfoFilter(logging.Filter):
        """Allows only INFO level log records."""

        def filter(self, record: logging.LogRecord) -> bool:
            return record.levelno == logging.INFO

    class ElseFilter(logging.Filter):
        """Allows all log records except INFO level."""

        def filter(self, record: logging.LogRecord) -> bool:
            return record.levelno != logging.INFO

    def __init__(
        self,
        debug: bool = False,
        custom_config: Optional[Dict] = None,
        ini_file: Optional[str] = None,
    ) -> None:
        """
        Initializes the LoggingConfigurator.

        Args:
            debug (bool): If True, uses the detailed formatter. If False, uses the simple formatter with no level name for INFO messages.
            custom_config (Optional[Dict]): Dictionary containing custom configurations for logging.
            ini_file (Optional[str]): Path to an INI file containing additional logging configurations.
        """
        self.debug = debug
        self.custom_config = custom_config or {}
        self.ini_file = ini_file
        self.configure()

    @staticmethod
    def get_logger(
        debug: bool = False,
        custom_config: Optional[Dict] = None,
        ini_file: Optional[str] = None,
    ) -> logging.Logger:
        """Retrieves the logger configured for 'covsonar'.

        Args:
            debug (bool): If True, uses the detailed formatter. If False, uses the simple formatter with no level name for INFO messages.
            custom_config (Optional[Dict]): Dictionary containing custom configurations for logging.
            ini_file (Optional[str]): Path to an INI file containing additional logging configurations.

        Returns:
            logging.Logger: The configured logger instance.
        """
        logger = logging.getLogger("covsonar")

        if not logger.hasHandlers():
            LoggingConfigurator(
                debug=debug, custom_config=custom_config, ini_file=ini_file
            )

        return logger

    def remove_logger_config(self):
        """
        Removes the logger config.
        """
        logger = logging.getLogger("covsonar")

        # Remove all handlers associated with the logger.
        for handler in logger.handlers[:]:
            handler.close()
            logger.removeHandler(handler)

    def set_debug_mode(self, debug: bool) -> None:
        """
        Enables or disables debug mode for logging.

        Args:
            debug (bool): If True, enables debug mode with detailed formatting. If False, disables debug mode with simple formatting.
        """
        self.debug = debug
        self.configure()  # Reconfigure the logging based on the new debug value

    def configure(self) -> None:
        """Configures the logging based on the debug flag, optional custom configurations, and optional INI file."""
        logger = logging.getLogger("covsonar")

        # Remove all existing handlers
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)

        # Read from INI file if provided
        if self.ini_file:
            config_parser = configparser.ConfigParser(interpolation=None)
            config_parser.read(self.ini_file)
            for section_items in config_parser.values():
                self.custom_config.update(section_items)

        else_handler = logging.StreamHandler(sys.stderr)
        else_handler.setLevel(logging.DEBUG)
        else_handler.addFilter(self.ElseFilter())

        info_handler = logging.StreamHandler(sys.stdout)
        info_handler.setLevel(logging.INFO)
        info_handler.addFilter(self.InfoFilter())

        if self.debug:
            logger.setLevel(logging.DEBUG)
            formatter_detailed = logging.Formatter(
                self.custom_config.get(
                    "detailed_format",
                    "%(asctime)s - %(name)s - %(filename)s:%(lineno)d - %(levelname)s: %(message)s",
                )
            )
            info_handler.setFormatter(formatter_detailed)
            else_handler.setFormatter(formatter_detailed)
        else:
            logger.setLevel(logging.INFO)
            formatter_info = logging.Formatter("%(message)s")
            formatter_else = logging.Formatter("%(levelname)s: %(message)s")

            info_handler.setFormatter(formatter_info)
            else_handler.setFormatter(formatter_else)

        logger.addHandler(info_handler)
        logger.addHandler(else_handler)
