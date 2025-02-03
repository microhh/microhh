import logging
from colorlog import ColoredFormatter

logger = logging.getLogger("puhhpy")
logger.setLevel(logging.DEBUG)

if not logger.handlers:
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.DEBUG)

    formatter = ColoredFormatter(
        "[%(asctime)s] [%(name)s] %(log_color)s[%(levelname)s] %(message)s'\033[0m",
        datefmt="%Y/%m/%d %H:%M:%S",
        log_colors={
            "DEBUG": "fg_244",
            "INFO": "",  # No color for INFO
            "WARNING": "fg_208",
            "ERROR": "red",
            "CRITICAL": "red",
        },
    )

    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    def critical_exception(message, *args, **kwargs):
         logger._log(logging.CRITICAL, message, args, **kwargs)
         raise RuntimeError(message)

    logger.critical = critical_exception